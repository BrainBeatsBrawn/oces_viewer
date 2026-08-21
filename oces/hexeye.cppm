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
            this->init (r_eye);
        }

        // @refeye: the reference OCES eye from which we are created.
        void init (const oces::eye& refeye)
        {
            // Find/set boundary for the hexgrid that matches the eye defined in the oces::reader.

            // Set up hex grid. Need meabn ommatidial neighbour distance to create the hexgrid, along with the approximate span.
            this->hg.init (refeye.d_mean, refeye.d_max * 1.5f);
        }
    };
}
