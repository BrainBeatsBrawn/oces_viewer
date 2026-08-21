// -*- C++ -*-
/*
 * Make a hexagonally arranged oces eye, with an sm::hexgrid used to define the eye layout.
 */
module;

#include <iostream>
#include <cstdint>

export module oces.hexeye;
export import oces.reader;

export import sm.hexgrid;
//export import sm.hexgrid_hdf; // probably
import sm.vvec;
import sm.mat;
import sm.bezcurvepath;
import sm.bezcurve;
import sm.centroid;

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

            // start with the z values
            sm::vvec<float> eye_z (refeye.eye_plane_coordinates.size(), 0.0f);
            for (std::uint32_t i = 0; i < refeye.eye_plane_coordinates.size(); ++i) {
                eye_z[i] = refeye.eye_plane_coordinates[i][2];
                coords2[i][0] = refeye.eye_plane_coordinates[i][0];
                coords2[i][1] = refeye.eye_plane_coordinates[i][1]; // already populated really
            }

            this->hg_z = this->hg.resample_data (eye_z, coords2, refeye.d_mean);
            std::cout << "hg_z: " << hg_z << std::endl;
        }
    };
}
