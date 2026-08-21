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

            // Find most central ommatidium and its neighbours
#if 0
            // central omm:
            refeye.eye_plane_coordinates[refeye.central_omm];
            // neigbours
            for (std::uint32_t k = 0; k < 6; ++k) {
                refeye.central_neighbours[k];
            }
#endif

            // Find/set boundary for the hexgrid that matches the eye defined in the oces::reader.

            // Set up hex grid. Need mean ommatidial neighbour distance to create the hexgrid, along with the approximate span.
            this->hg.init (refeye.d_mean, refeye.d_max * 2.0f);

            auto hiter = this->hg.hexen.begin();
            std::cout << "Hexgrid centre is at 0? " << hiter->output_xy() << std::endl;

            // Now find the angles...
            float ang = refeye.central_neighbour_angles.min();
            sm::mat<float, 4> tfm (sm::quaternion<float>(sm::vec<float>::uz(), -ang));
            sm::vec<float> c = refeye.eye_plane_coordinates[refeye.central_omm];
            c[2] = 0;
            tfm.pretranslate (-c);
            this->hg.transform (tfm);

            this->hg.set_circular_boundary (refeye.d_max * 1.2f / 2.0f);
        }
    };
}
