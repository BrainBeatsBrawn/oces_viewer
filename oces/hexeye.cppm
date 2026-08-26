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
        // The hexy eye we create
        oces::eye eye;
        // The flat hexgrid
        sm::hexgrid<F> hg;
        // The 'z' positions of the hexgrid eye
        sm::vvec<F> hg_z;
        // The boundary around the outer-most ommatidia. The eye outline.
        sm::bezcurvepath<F> eye_outline;

        hexeye(){}
        // @refeye: the reference OCES eye from which we are created.
        hexeye (const oces::eye& refeye) { this->init (refeye); }

        // Average height from any point within one or two hexes of the current hex
        sm::vvec<F> use_average_nearest_z (const sm::vvec<F>& eye_z,
                                           const sm::vvec<sm::vec<F, 2>>& coords2)
        {
            sm::vvec<F> hex_z (this->hg.num(), 0.0f);

            const F d_thresh = this->hg.d * F{1};

            // For each hex, find all the close real data and copy the z values/average.
            for (std::uint32_t xi = 0u; xi < this->hg.num(); ++xi) {

                const sm::vec<F, 2> hp = { this->hg.d_x[xi], this->hg.d_y[xi] };

                sm::vvec<std::uint32_t> nearby = {};

                for (std::uint32_t i = 0u; i < coords2.size(); ++i) {
                    F _d = (hp - coords2[i]).length(); // distance between coord and our hex
                    if (_d < d_thresh) { nearby.push_back (i); }
                }
                // Found nearby, now use them
                F z_sum = F{0};
                for (std::uint32_t i = 0u; i < nearby.size(); ++i) {
                    z_sum += eye_z[nearby[i]];
                }
                hex_z[xi] = z_sum / nearby.size();
            }
            return hex_z;
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
            F ang = refeye.central_neighbour_angles.min();
            sm::mat<F, 4> tfm (sm::quaternion<F>(sm::vec<F>::uz(), ang));
            sm::vec<F> c = refeye.eye_plane_coordinates[refeye.central_omm];
            c[2] = 0;
            tfm.pretranslate (c);
            this->hg.transform (tfm);

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
                this->eye_outline.add_curve (crv);
            }
            // Add the closing path
            s = bnd2[bnd2.size() - 1];
            e = bnd2[0];
            sm::bezcurve<F, 3> crv (s, e, e, s); // third order, but make it a straight section from bnd2.
            this->eye_outline.add_curve (crv);


            // When we set boundary, we haven't closed the loop.
            this->hg.set_boundary (eye_outline);

            // Now resample eye_plane_coordinates[][2], refeye.directions and refeye.acceptance_angles onto hex containers
            sm::vvec<sm::vec<F>> all_coordinates (refeye.eye_plane_coordinates.size());
            for (std::uint32_t i = 0; i < refeye.eye_plane_coordinates.size(); ++i) {
                all_coordinates[i] = refeye.eye_plane_coordinates[i].template as<F>();
            }

            // start with the z values
            sm::vvec<F> eye_z (all_coordinates.size(), 0.0f);
            coords2.resize (all_coordinates.size());
            for (std::uint32_t i = 0; i < all_coordinates.size(); ++i) {
                const sm::vec<F> crd = refeye.eye_plane_coordinates[i].template as<F>();
                coords2[i][0] = crd[0];
                coords2[i][1] = crd[1];
                eye_z[i]      = crd[2];
            }

            // this->hg_z = this->hg.resample_data (eye_z, coords2, 2 * refeye.d_mean);
            this->hg_z = this->use_average_nearest_z (eye_z, coords2);
        }
    };
}
