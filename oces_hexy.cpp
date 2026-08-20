/*
 * Read OCES file then display with mathplot::compoundray::EyeVisual (no 2D projection here)
 *
 * THEN create an equivalent HexGrid eye, in which neighbour relationships between ommatidia are
 * known, and which has the right outline and ommatidia directions.
 */

#include <memory>
#include <iostream>
#include <string>
#define ARGS_NOEXCEPT 1
#include <args/args.hxx> // github.com/Taywee/args

import oces.reader;
import mplot.visual;
import mplot.spherevisual;
import mplot.colourbarvisual;
import mplot.compoundray.eyevisual;

int main (int argc, char** argv)
{
    args::ArgumentParser ap ("OCES viewer", "Have a nice day.");
    args::ValueFlag<std::string> a_fname  (ap, "filepath", "path/to/oces_file.gltf",       {'f'});
    args::Flag a_fov (ap, "fov", "Show field of view with acceptance angle cones", {'v', "fov"});
    args::Flag a_hidehead (ap, "hidehead", "Hide the head, even if it was read from OCES file", {'i', "hidehead"});
    args::Flag a_showrays (ap, "showrays", "Show the 2D ommatidia projection rays", {'y', "showrays"});
    ap.ParseCLI (argc, argv);
    std::string filename = "";
    if (a_fname) {
        filename = args::get (a_fname);
    } else {
        std::cerr << ap;
        return -1;
    }

    // Read
    oces::reader oces_reader (filename);
    if (oces_reader.read_success == false) { return -1; }
    // Now view
    auto v = mplot::Visual<>(1024, 768, "mplot::compoundray::EyeVisual");
    v.lightingEffects (true);
    // We read the information from the eye file into a vector of Ommatidium objects.  Ommatidium is
    // defined in "cameras/CompoundEyeDataTypes.h" in compound ray, mplot::Ommatidium is a
    // mplot/Seb's maths style equivalent. It contains 2 3D float vectors and two scalar floating point
    // values.
    auto ommatidia = std::make_unique<std::vector<mplot::compoundray::Ommatidium>>();
    std::vector<std::array<float, 3>> ommatidiaColours;
    // Copy data into the ommatidia data structure
    ommatidia->resize (oces_reader.eye.position.size());
    std::cerr << "Copying " << oces_reader.eye.position.size() << " ommatidia\n";
    for (size_t i = 0; i < oces_reader.eye.position.size(); ++i) {
        (*ommatidia)[i].relativePosition = oces_reader.eye.position[i];
        (*ommatidia)[i].relativeDirection = oces_reader.eye.orientation[i];
        (*ommatidia)[i].focalPointOffset = oces_reader.eye.focal_offset[i];
        (*ommatidia)[i].acceptanceAngleRadians = oces_reader.eye.acceptance_angle[i];
    }
    // Make some dummy data to demo the eye
    sm::vvec<float> ommatidiaData;
    ommatidiaData.linspace (0, 1, ommatidia->size());
    // Colour it with a colour map
    mplot::ColourMap cm (mplot::ColourMapType::Plasma);
    ommatidiaColours.resize (ommatidia->size());
    for (size_t i = 0; i < ommatidia->size(); ++i) { ommatidiaColours[i] = cm.convert (ommatidiaData[i]); }

    mplot::meshgroup* head_mesh_ptr = nullptr;
    if (!a_hidehead) { head_mesh_ptr = reinterpret_cast<mplot::meshgroup*>(&oces_reader.eye.head_mesh); }

    oces_reader.eye.head_mesh.single_colour = {0.345f, 0.122f, 0.082f};
    auto eyevm = std::make_unique<mplot::compoundray::EyeVisual<>> (sm::vec<>{}, &ommatidiaColours, ommatidia.get(), head_mesh_ptr);
    eyevm->set_parent (v.get_id());
    eyevm->name = "CompoundRay Eye";
    eyevm->show_cones = true;
    eyevm->show_rays = (a_showrays ? true : false);
    eyevm->show_cones = false;
    eyevm->show_fov = (a_fov ? true : false);
    eyevm->pre_set_cone_length (0.005f);
    eyevm->finalize();
    auto ep = v.addVisualModel (eyevm);
    ep->scaleViewMatrix (1000.0f);

    // Now find that HexGrid.

    // What's the mean position of the ommatidia? What's the 'normal' of the eye? We'll start with a
    // flat HexGrid whose plane is defined by this normal. Best fit ellipsoid?

    v.keepOpen();

    return 0;
}
