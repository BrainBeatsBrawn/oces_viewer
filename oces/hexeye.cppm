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
import sm.bezcurvepath;
import sm.bezcurve;
import sm.centroid;
import sm.algo;

export namespace oces
{
    struct hexeye
    {
        oces::eye eye;    // The hexy eye we create
        sm::hexgrid hg;   // The flat hexgrid
        sm::vvec<float> hg_z; // The 'z' positions of the hexgrid eye

        hexeye(){}

        // @refeye: the reference OCES eye from which we are created.
        hexeye (const oces::eye& refeye)
        {
            this->init (refeye);
        }

        void find_extra_coordinates (const oces::eye& refeye, [[maybe_unused]] sm::vvec<sm::vec<float>>& extra_coordinates)
        {
            const std::uint32_t omm_per_eye = refeye.eye_plane_coordinates.size();

            // Map by angle
            std::map<std::uint32_t, std::uint32_t> idx_by_angle; // key: angle index, value: coordinate index
            std::map<std::uint32_t, float> dist_by_angle;
            std::uint32_t n_angles = 128;
            for (std::uint32_t i = 0; i < omm_per_eye; ++i) {
                sm::vec<float, 2> c2 = refeye.eye_plane_coordinates[i].less_one_dim();
                // Use c2.angle to bin this into n_angles bins. If it is the furthest for that
                // angle, then record in idx_by_angle.
                float d = c2.length();
                float a = c2.angle();
                sm::algo::zero_to_twopi (a);
                std::uint32_t aidx = static_cast<std::uint32_t>(n_angles * a / sm::mathconst<float>::two_pi);
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
                sm::vec<float, 3> nearest_3 = {};
                nearest_3.set_from (std::numeric_limits<float>::max());

                sm::vec<std::uint32_t, 3> nearest_3_cidx;

                const auto& pi = refeye.eye_plane_coordinates[cidx];

                for (std::uint32_t j = 0; j < omm_per_eye; ++j) {
                    if (j == cidx) { continue; }
                    const auto& pj = refeye.eye_plane_coordinates[j];
                    const float dij = (pi - pj).length();
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

                std::cout << "Have pairs "
                          << cidx << "-" << nearest_3_cidx[0] << ", "
                          << cidx << "-" << nearest_3_cidx[1] << ", "
                          << cidx << "-" << nearest_3_cidx[2] << "\n";
            }
        }

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
            float ang = refeye.central_neighbour_angles.min();
            sm::mat<float, 4> tfm (sm::quaternion<float>(sm::vec<float>::uz(), -ang));
            sm::vec<float> c = refeye.eye_plane_coordinates[refeye.central_omm];
            c[2] = 0;
            tfm.pretranslate (-c);
            this->hg.transform (tfm);

            // Need 2D points for graham scan to find the convex hull of eye_plane_coordinates
            sm::vvec<sm::vec<float, 2>> coords2 (refeye.get_omm_per_eye());
            for (std::uint32_t i = 0; i < refeye.get_omm_per_eye(); ++i) {
                coords2[i] = refeye.eye_plane_coordinates[i].less_one_dim();
                coords2[i][0] *= -1; // Not quite sure why I have to invert x.
            }
            sm::vvec<sm::vec<float, 2>> bnd2 = sm::geometry::graham_scan (coords2);
            // From bnd2 make up a set of bezcoords to form a bezcurvepath
            sm::bezcurvepath<float> bcp;
            for (std::uint32_t i = 1; i < bnd2.size(); ++i) {
                sm::vec<float, 2> s = bnd2[i-1];
                sm::vec<float, 2> e = bnd2[i];
                sm::bezcurve<float, 3> crv (s, e, e, s); // third order, but make it a straight section from bnd2.
                bcp.add_curve (crv);
            }
            this->hg.set_boundary (bcp);

            // Now resample eye_plane_coordinates[][2], refeye.directions and refeye.acceptance_angles onto hex containers

            // To get a good resampling, it helps to extend the data in eye_plane_coordinates so
            // that the data extends beyond the boundary. I'll extend with a linear extension. Find
            // points around the outside of the eye. For each found point, find a nearby point, and
            // draw a line between the two beyond the boundary.

            sm::vvec<sm::vec<float>> extra_coordinates;
            this->find_extra_coordinates (refeye, extra_coordinates); // may also pass in dirns, acceptance angles

            sm::vvec<sm::vec<float>> all_coordinates = refeye.eye_plane_coordinates;
            all_coordinates.append (extra_coordinates);

            // start with the z values
            sm::vvec<float> eye_z (all_coordinates.size(), 0.0f);
            coords2.resize (all_coordinates.size());
            for (std::uint32_t i = 0; i < all_coordinates.size(); ++i) {
                coords2[i][0] = all_coordinates[i][0];
                coords2[i][1] = all_coordinates[i][1];
                eye_z[i]      = all_coordinates[i][2];
            }

            this->hg_z = this->hg.resample_data (eye_z, coords2, refeye.d_mean);
            std::cout << "hg_z: " << hg_z << std::endl;
        }
    };
}
