// -*- C++ -*-
/*
 * Make a hexagonally arranged oces eye, with an sm::hexgrid used to define the eye layout.
 */
module;

#include <iostream>
#include <cstdint>
#include <map>

export module oces.hexeye;
export import oces.reader;

export import sm.hexgrid;
//export import sm.hexgrid_hdf; // probably
import sm.vvec;
import sm.mat;
export import sm.bezcurvepath;
import sm.bezcurve;
import sm.centroid;
import sm.algo;

export namespace oces
{
    template<typename F=float>
    struct hexeye
    {
        oces::eye eye;    // The hexy eye we create
        sm::hexgrid<F> hg;   // The flat hexgrid
        sm::vvec<F> hg_z; // The 'z' positions of the hexgrid eye

        hexeye(){}

        // @refeye: the reference OCES eye from which we are created.
        hexeye (const oces::eye& refeye)
        {
            this->init (refeye);
        }

        void find_extra_coordinates (const oces::eye& refeye, [[maybe_unused]] sm::vvec<sm::vec<F>>& extra_coordinates)
        {
            const std::uint32_t omm_per_eye = refeye.eye_plane_coordinates.size();

            // Map by angle
            std::map<std::uint32_t, std::uint32_t> idx_by_angle; // key: angle index, value: coordinate index
            std::map<std::uint32_t, F> dist_by_angle;
            std::uint32_t n_angles = 128;
            for (std::uint32_t i = 0; i < omm_per_eye; ++i) {
                sm::vec<F, 2> c2 = refeye.eye_plane_coordinates[i].less_one_dim();
                // Use c2.angle to bin this into n_angles bins. If it is the furthest for that
                // angle, then record in idx_by_angle.
                F d = c2.length();
                F a = c2.angle();
                sm::algo::zero_to_twopi (a);
                std::uint32_t aidx = static_cast<std::uint32_t>(n_angles * a / sm::mathconst<F>::two_pi);
                if (dist_by_angle.contains (aidx)) {
                    if (d > dist_by_angle[aidx]) {
                        dist_by_angle[aidx] = d;
                        idx_by_angle[aidx] = i;
                    }
                } else {
                    dist_by_angle[aidx] = d;
                    idx_by_angle[aidx] = i;
                }
            }

            std::cout << "idx_by_angle size:"  << idx_by_angle.size() << std::endl;

            for (auto const& [aidx, cidx] : idx_by_angle) {

                // Now, for each coordinate indexed in idx_by_angle, find its 2 or 3 nearest neighbours, and draw new points
                sm::vec<F, 3> nearest_3 = {};
                nearest_3.set_from (std::numeric_limits<F>::max());

                sm::vec<std::uint32_t, 3> nearest_3_cidx;

                const auto& pi = refeye.eye_plane_coordinates[cidx];

                for (std::uint32_t j = 0; j < omm_per_eye; ++j) {
                    if (j == cidx) { continue; }
                    const auto& pj = refeye.eye_plane_coordinates[j];
                    const F dij = (pi - pj).length();
                    // If dij is less than any, replace the biggest
                    auto kbig = nearest_3.argmax();
                    for (std::uint32_t k = 0; k < 3; ++k) {
                        if (dij <= nearest_3[k]) {
                            nearest_3[kbig] = dij;
                            nearest_3_cidx[kbig] = j;
                            break;
                        }
                    }
                }

                //std::cout << "Have pairs "
                //          << cidx << "-" << nearest_3_cidx[0] << ", "
                //          << cidx << "-" << nearest_3_cidx[1] << ", "
                //          << cidx << "-" << nearest_3_cidx[2] << "\n";

                for (std::uint32_t k = 0; k < 3; ++k) {
                    // pi and nearest_3_cidx[k]
                    sm::vec<F> c = refeye.eye_plane_coordinates[nearest_3_cidx[k]];
                    // Find order
                    sm::vec<F> pi_c = c - pi;
                    sm::vec<F> c_pi = -pi_c;
                    sm::vec<F> new_p = {};
                    sm::vec<F> new_p2 = {};
                    if (pi_c.dot (pi) > 0.0f) {
                        // pi_c is outward
                        new_p = pi + pi_c * 2.0f;
                        new_p2 = pi + pi_c * 4.0f;
                    } else {
                        // c_pi outward
                        new_p = pi + c_pi * 2.0f;
                        new_p2 = pi + c_pi * 4.0f;
                    }
                    extra_coordinates.push_back (new_p);
                    extra_coordinates.push_back (new_p2);
                    // Plus work out dirn, acceptance angle etc
                }
            }
        }

        // Just find the closest height from eye_z and place that in the hexgrid.
        sm::vvec<F> use_nearest_z (const sm::vvec<F>& eye_z,
                                   const sm::vvec<sm::vec<F, 2>>& coords2)
        {
            sm::vvec<F> hex_z (this->hg.num(), 0.0f);

            // For each hex, find the closest real datum and copy the z value.
            for (std::uint32_t xi = 0u; xi < this->hg.num(); ++xi) {
                const sm::vec<F, 2> hp = { this->hg.d_x[xi], this->hg.d_y[xi] };
                std::uint32_t closest = 0u;
                F d_closest = std::numeric_limits<F>::max();
                for (std::uint32_t i = 0u; i < coords2.size(); ++i) {
                    F _d = (hp - coords2[i]).length(); // distance between coord and our hex
                    if (_d < d_closest) {
                        d_closest = _d;
                        closest = i;
                    }
                }
                // Found closest, now use it
                hex_z[xi] = eye_z[closest];
            }
            return hex_z;
        }

        sm::vvec<sm::vec<F>> extra_coordinates;

        sm::bezcurvepath<F> bcp;

        sm::interval<F> eye_z_range;

        // @refeye: the reference OCES eye from which we are created.
        void init (const oces::eye& refeye)
        {
            if (refeye.ready == false) {
                std::cerr << "Can't use that reference eye as it is not ready. Returning.\n";
                return;
            }

            // Set up hex grid. Need mean ommatidial neighbour distance to create the hexgrid, along with the approximate span.
            this->hg.init (refeye.d_mean, refeye.d_max * 2.0f);

            // Now use offset and angle to make a transform for the hexgrid
#if 0
            F ang = refeye.central_neighbour_angles.min();
            sm::mat<F, 4> tfm (sm::quaternion<F>(sm::vec<F>::uz(), -ang));
            sm::vec<F> c = refeye.eye_plane_coordinates[refeye.central_omm];
            c[2] = 0;
            tfm.pretranslate (-c);
            this->hg.transform (tfm);
#endif
            // Need 2D points for graham scan to find the convex hull of eye_plane_coordinates
            sm::vvec<sm::vec<F, 2>> coords2 (refeye.get_omm_per_eye());
            for (std::uint32_t i = 0; i < refeye.get_omm_per_eye(); ++i) {
                coords2[i] = refeye.eye_plane_coordinates[i].less_one_dim().template as<F>();
                //coords2[i][1] *= -1; // invert y. I also invert y when I plot with BezCurvePathVisual. Figure out why.
            }
            sm::vvec<sm::vec<F, 2>> bnd2 = sm::geometry::graham_scan (coords2);
            sm::vec<F, 2> s = {};
            sm::vec<F, 2> e = {};
            // From bnd2 make up a set of bezcoords to form a bezcurvepath
            for (std::uint32_t i = 1; i < bnd2.size(); ++i) {
                s = bnd2[i-1];
                e = bnd2[i];
                sm::bezcurve<F, 3> crv (s, e, e, s); // third order, but make it a straight section from bnd2.
                this->bcp.add_curve (crv);
            }
            // Add the closing path
            s = bnd2[bnd2.size() - 1];
            e = bnd2[0];
            sm::bezcurve<F, 3> crv (s, e, e, s); // third order, but make it a straight section from bnd2.
            this->bcp.add_curve (crv);


            // When we set boundary, we haven't closed the loop.
            this->hg.set_boundary (bcp);

            // Now resample eye_plane_coordinates[][2], refeye.directions and refeye.acceptance_angles onto hex containers

            // To get a good resampling, it may help to extend the data in eye_plane_coordinates so
            // that the data extends beyond the boundary. I'll extend with a linear extension. Find
            // points around the outside of the eye. For each found point, find a nearby point, and
            // draw a line between the two beyond the boundary.
            // this->find_extra_coordinates (refeye, this->extra_coordinates); // may also pass in dirns, acceptance angles

            sm::vvec<sm::vec<F>> all_coordinates (refeye.eye_plane_coordinates.size());
            for (std::uint32_t i = 0; i < refeye.eye_plane_coordinates.size(); ++i) {
                all_coordinates[i] = refeye.eye_plane_coordinates[i].template as<F>();
            }
            std::cout << "eye_plane_coordinates size " << refeye.eye_plane_coordinates.size() << std::endl;
            all_coordinates.append (this->extra_coordinates);
            std::cout << "extra_coordinates size " << this->extra_coordinates.size() << std::endl;

            sm::vvec<F> eye_z0 (refeye.eye_plane_coordinates.size(), 0.0f);
            for (std::uint32_t i = 0; i < refeye.eye_plane_coordinates.size(); ++i) {
                eye_z0[i] = refeye.eye_plane_coordinates[i][2];
            }
            this->eye_z_range = eye_z0.range();

            // start with the z values
            sm::vvec<F> eye_z (all_coordinates.size(), 0.0f);
            coords2.resize (all_coordinates.size());
            for (std::uint32_t i = 0; i < all_coordinates.size(); ++i) {
                coords2[i][0] = all_coordinates[i][0];
                coords2[i][1] = all_coordinates[i][1];
                eye_z[i]      = all_coordinates[i][2];
            }

#if 1
            this->hg_z = this->hg.resample_data (eye_z, coords2, 2 * refeye.d_mean);
#else
            this->hg_z = this->use_nearest_z (eye_z, coords2);
#endif
            //std::cout << "hg_z: " << hg_z << std::endl;
        }
    };
}
