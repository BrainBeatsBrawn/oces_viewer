/*
 * Read OCES file then display with mathplot::compoundray::EyeVisual
 */

#include <memory>
#include <format>
#include <iostream>
#include <string>
#define ARGS_NOEXCEPT 1
#include <args/args.hxx> // github.com/Taywee/args

import oces.reader;
import mplot.visual;
import mplot.spherevisual;
import mplot.colourbarvisual;
import mplot.coordarrows;
import mplot.quivervisual;
import mplot.compoundray.eyevisual;
import sm.flags;

// Enumerate our view options
enum class viewopts : std::uint32_t
{
    show_fov,
    show_proj_fov,
    show_sphere,
    show_twod,
    show_rays,
    show_head,
    show_origin
};

// Extend mplot::Visual for key-commands
struct OcesVisual final : public mplot::Visual<>
{
    // Boilerplate constructor calls Visual constructor
    OcesVisual(int width, int height, const std::string & title) : mplot::Visual<>(width, height, title) {}
    // Boolean flags for what we want to have visible
    sm::flags<viewopts> view_options;
protected:
    // Add our key event handling in this extension of mplot::Visual::key_callback_extra
    void key_callback_extra ([[maybe_unused]] int key, [[maybe_unused]] int scancode,
                             [[maybe_unused]] int action, [[maybe_unused]] int mods) override
    {
        if (key == mplot::key::v && action == mplot::keyaction::press) { this->view_options.flip (viewopts::show_fov); }
        if (key == mplot::key::f && action == mplot::keyaction::press) { this->view_options.flip (viewopts::show_proj_fov); }
        if (key == mplot::key::s && action == mplot::keyaction::press) { this->view_options.flip (viewopts::show_sphere); }
        if (key == mplot::key::p && action == mplot::keyaction::press) { this->view_options.flip (viewopts::show_twod); }
        if (key == mplot::key::y && action == mplot::keyaction::press) { this->view_options.flip (viewopts::show_rays); }
        if (key == mplot::key::i && action == mplot::keyaction::press) { this->view_options.flip (viewopts::show_head); }
        if (key == mplot::key::z && action == mplot::keyaction::press) { this->view_options.flip (viewopts::show_origin); }
        // Add app-specific help output:
        if (key == mplot::key::h && action == mplot::keyaction::press) {
            std::cout << "OCES Viewer key commands:\n";
            std::cout << "v: Toggle FOV\n";
            std::cout << "f: Toggle projected horz/vert FOV\n";
            std::cout << "s: Toggle Sphere\n";
            std::cout << "p: Toggle 2D projection of eye for visualizations\n";
            std::cout << "y: Toggle Rays\n";
            std::cout << "i: Toggle Head\n";
            std::cout << "z: Toggle OCES origin\n";
        }
    }
};

// Helper to make (and remake) the eye model
mplot::compoundray::EyeVisual<>* make_eye_model (OcesVisual& v, oces::reader& oces_reader,
                                                 std::vector<mplot::compoundray::Ommatidium>* ommatidia,
                                                 std::vector<std::array<float, 3>>* ommatidiaColours,
                                                 const std::string& projstr,
                                                 const float psrad, sm::vec<float> pscentre, const float psr,
                                                 const sm::vec<float>& psrax, const sm::vec<float>& twod_shift,
                                                 mplot::compoundray::EyeVisual<>* ep)
{
    if (ep != nullptr) { v.removeVisualModel (ep); }

    mplot::meshgroup* head_mesh_ptr = nullptr;
    if (v.view_options.test (viewopts::show_head)) {
        head_mesh_ptr = reinterpret_cast<mplot::meshgroup*>(&oces_reader.head_mesh);
    }

    // Set a fixed colour for head mesh
    oces_reader.head_mesh.single_colour = {0.345f, 0.122f, 0.082f};

    auto eyevm = std::make_unique<mplot::compoundray::EyeVisual<>> (sm::vec<>{}, ommatidiaColours, ommatidia, head_mesh_ptr);
    eyevm->set_parent (v.get_id());
    eyevm->name = "CompoundRay Eye";
    eyevm->show_cones = true;

    [[maybe_unused]] auto ptype = mplot::compoundray::EyeVisual<>::projection_type::equirectangular; // mercator, equirectangular or cassini
    if (projstr.find ("merc") != std::string::npos) {
        ptype = mplot::compoundray::EyeVisual<>::projection_type::mercator;
    } else if (projstr.find ("cass") != std::string::npos) {
        ptype = mplot::compoundray::EyeVisual<>::projection_type::cassini;
    }

    sm::mat<float, 4> twod_tr;
    twod_tr.translate (twod_shift);

    if (psrad > 0.0f && v.view_options.test (viewopts::show_twod)) {
        // To avoid 2D, don't add spherical projections
        std::cout << "Rotation about axis " << psrax << " by amount " << psr << " radians\n";
        sm::quaternion<float> psrotn (psrax, psr);
        eyevm->add_spherical_projection (ptype, twod_tr, pscentre, psrad, psrotn, 0, oces_reader.position.size() / 2);
        if (oces_reader.mirrors.empty() == false) {
            pscentre = (oces_reader.mirrors[0] * pscentre).less_one_dim();
            std::cout << "New centre: " << pscentre << std::endl;

            sm::vec<> twod_shift_left = twod_shift;
            twod_shift_left[0] *= -1.0f;
            twod_tr.set_identity();
            twod_tr.translate (twod_shift_left);
            eyevm->add_spherical_projection (ptype, twod_tr, pscentre, psrad, psrotn.invert(),
                                             oces_reader.position.size() / 2, oces_reader.position.size());
        }
    }
    eyevm->show_sphere = v.view_options.test (viewopts::show_sphere);
    eyevm->show_rays = v.view_options.test (viewopts::show_rays);
    eyevm->show_cones = false;
    eyevm->show_fov = v.view_options.test (viewopts::show_fov);
    eyevm->pre_set_cone_length (0.005f);
    eyevm->finalize();

    mplot::compoundray::EyeVisual<>* _ep = v.addVisualModel (eyevm);
    _ep->scaleViewMatrix (1000.0f); // Tiny ant eyes are scaled by a big factor to be in more useable model units
    return _ep;
}

int main (int argc, char** argv)
{
    args::ArgumentParser ap ("OCES viewer", "Have a nice day.");

    args::ValueFlag<std::string> a_fname  (ap, "filepath", "path/to/oces_file.gltf",       {'f'});
    args::ValueFlag<float>       a_psrad  (ap, "radius",   "The projection sphere radius (numeric value)", {'r'});
    args::ValueFlag<std::string> a_centre (ap, "centre",   "The projection sphere centre (comma separated coordinates)", {'c'});
    args::ValueFlag<std::string> a_psrax (ap, "psrax",   "The projection sphere rotation axis (comma separated coordinates)", {'x'});
    args::ValueFlag<std::string> a_2dshift (ap, "twodshift",   "A translation of the 2d right eye (comma separated coordinates, mirrored for left eye)", {'d'});
    args::ValueFlag<float>       a_psr  (ap, "psr",   "The projection sphere rotation radians (numeric value)", {'o'});
    args::ValueFlag<std::string> a_proj (ap, "proj",   "The projection type (equirectangular, mercator or cassini)", {'p'});
    args::Flag a_fov (ap, "fov", "Show field of view with acceptance angle cones", {'v', "fov"});
    args::Flag a_verbose (ap, "verbose", "Print ommatidial positions and orientations to stdout", {'V', "verbose"});
    args::Flag a_hidemirror (ap, "hidemirror", "Hide the mirrored eye", {'H', "hidemirror"});
    args::Flag a_hidehead (ap, "hidehead", "Hide the head, even if it was read from OCES file", {'i', "hidehead"});
    args::Flag a_showsphere (ap, "showsphere", "Show the 2D projection sphere", {'s', "showsphere"});
    args::Flag a_showrays (ap, "showrays", "Show the 2D ommatidia projection rays", {'y', "showrays"});
    args::Flag a_showorigin (ap, "showorigin", "Show coordinate arrows in-scene at the OCES origin", {'z', "showorigin"});

    ap.ParseCLI (argc, argv);

    std::string filename = "";
    float psrad = 0.1f;
    sm::vec<float> pscentre = { 0, 0, 0 };
    float psr = 0.0f;
    sm::vec<float> psrax = { 0, 1, 0 };
    sm::vec<float> twod_shift = { -0.0005, 0.002, 0 }; // A shift of the twod representation of the
                                                       // right eye (the purple one in the example
                                                       // velox-head file)

    if (a_fname) {
        filename = args::get (a_fname);
    } else {
        std::cerr << ap;
        return -1;
    }

    if (a_psrad) {
        psrad = args::get (a_psrad);
        std::cerr << "User-supplied projection sphere radius: " << psrad << std::endl;
    }

    if (a_centre) {
        pscentre.set_from (args::get (a_centre));
        std::cerr << "User-supplied projection sphere centre: " << pscentre << std::endl;
    }

    if (a_psrax) {
        psrax.set_from (args::get (a_psrax));
        std::cerr << "User-supplied projection rotation axis: " << psrax << std::endl;
    }

    if (a_2dshift) {
        twod_shift.set_from (args::get (a_2dshift));
        std::cerr << "User-supplied 2D eye shift: " << twod_shift << std::endl;
    }

    if (a_psr) {
        psr = args::get (a_psr);
        std::cerr << "User-supplied projection rotation radians: " << psr << std::endl;
    }

    std::string projstr = "equirectangular";
    if (a_proj) {
        projstr = args::get (a_proj);
        std::cerr << "User-supplied projection type: " << projstr << std::endl;
    }

    // Read
    oces::reader oces_reader (filename, (a_hidemirror ? true : false));
    if (oces_reader.read_success == false) { return -1; }

    // Now view
    auto v = OcesVisual(1024, 768, "mplot::compoundray::EyeVisual");
    v.lightingEffects (true);

    // Copy cmd line options into v
    v.view_options.set (viewopts::show_fov, (a_fov ? true : false));
    v.view_options.set (viewopts::show_proj_fov, (a_fov ? true : false)); // use a_fov here, too
    v.view_options.set (viewopts::show_rays, (a_showrays ? true : false));
    v.view_options.set (viewopts::show_sphere, (a_showsphere ? true : false));
    v.view_options.set (viewopts::show_head, (a_hidehead ? false : true));
    v.view_options.set (viewopts::show_twod, (a_psrad ? true : false));
    v.view_options.set (viewopts::show_origin, (a_showorigin ? true : false));

    // We read the information from the eye file into a vector of Ommatidium objects.  Ommatidium is
    // defined in "cameras/CompoundEyeDataTypes.h" in compound ray, mplot::Ommatidium is a
    // mplot/Seb's maths style equivalent. It contains 2 3D float vectors and two scalar floating point
    // values.
    auto ommatidia = std::make_unique<std::vector<mplot::compoundray::Ommatidium>>();
    std::vector<std::array<float, 3>> ommatidiaColours;

    // Copy data into the ommatidia data structure
    ommatidia->resize (oces_reader.position.size());
    std::cerr << "Copying " << oces_reader.position.size() << " ommatidia\n";
    for (size_t i = 0; i < oces_reader.position.size(); ++i) {
        (*ommatidia)[i].relativePosition = oces_reader.position[i];
        (*ommatidia)[i].relativeDirection = oces_reader.orientation[i];
        (*ommatidia)[i].focalPointOffset = oces_reader.focal_offset[i];
        (*ommatidia)[i].acceptanceAngleRadians = oces_reader.acceptance_angle[i];
    }
    if (a_verbose) {
        for (size_t i = 0; i < oces_reader.position.size(); ++i) {
            std::cout << "ommatidium[" << i << "]: posn = " <<  oces_reader.position[i]
                      << ",  orientn = " << oces_reader.orientation[i]
                      << ",  fo = " << oces_reader.focal_offset[i]
                      << ",  aa = " << oces_reader.acceptance_angle[i]
                      << std::endl;
        }
    }

    // Make some dummy data to demo the eye
    sm::vvec<float> ommatidiaData;
    ommatidiaData.linspace (0, 1, ommatidia->size());
    // Colour it with a colour map
    mplot::ColourMap cm (mplot::ColourMapType::Plasma);
    ommatidiaColours.resize (ommatidia->size());
    for (size_t i = 0; i < ommatidia->size(); ++i) {
        ommatidiaColours[i] = cm.convert (ommatidiaData[i]);
    }

    // Place a colour bar for the ommtidia index
    auto cbv = std::make_unique<mplot::ColourBarVisual<float>>(sm::vec<>{0.5f,-0.6f,0});
    cbv->set_parent (v.get_id());
    cbv->orientation = mplot::colourbar_orientation::horizontal;
    cbv->tickside = mplot::colourbar_tickside::right_or_below;
    cbv->cm = cm;
    cbv->width = 0.08f;
    cbv->length = 0.5f;
    cbv->tf.fontsize = 0.04f;
    cbv->scale.compute_scaling (0.0f, static_cast<float>(ommatidia->size() - 1));
    cbv->addLabel ("Index 0-" + std::to_string (ommatidia->size() - 1),
                   sm::vec<>{0, 0.15f, 0}, mplot::TextFeatures(0.05f));
    cbv->finalize();
    auto cbvp = v.addVisualModel (cbv);
    cbvp->setHide (v.view_options.test (viewopts::show_fov) == true);

    // Coordinate arrows for OCES origin
    auto origin_ca = std::make_unique<mplot::CoordArrows<>> (sm::vec<>{});
    origin_ca->set_parent (v.get_id());
    origin_ca->em = 0.0f;
    float len = 2.0f;
    origin_ca->lengths = { len, len, len };
    origin_ca->thickness = 1.0f;
    origin_ca->endsphere_size = 1.2f;
    origin_ca->finalize();
    auto oces_origin = v.addVisualModel (origin_ca);
    oces_origin->setHide (v.view_options.test (viewopts::show_origin) == false);

    // Quivers for projected angles
    auto qvh = std::make_unique<mplot::QuiverVisual<float>>(&oces_reader.h_plane_position, sm::vec<>{-7, 0, 0},
                                                            &oces_reader.h_plane_orientation,
                                                            mplot::ColourMapType::MonochromeGreen);
    qvh->set_parent (v.get_id());
    qvh->quiver_length_gain = 4.0f;
    qvh->quiver_thickness_gain = 0.005f;
    qvh->addLabel (std::format ("Horz FOV: {:.2f}{}",
                                (oces_reader.horz_fov * sm::mathconst<float>::rad2deg),
                                mplot::unicode::toUtf8(mplot::unicode::degreesign)),
                   sm::vec<float>{-4,0,4}, mplot::TextFeatures(0.12f));
    qvh->finalize();
    auto qvhp_h = v.addVisualModel (qvh);
    qvhp_h->setHide (v.view_options.test (viewopts::show_proj_fov) == false);

    qvh = std::make_unique<mplot::QuiverVisual<float>>(&oces_reader.v_plane_position, sm::vec<>{-7, 0, 0},
                                                       &oces_reader.v_plane_orientation,
                                                       mplot::ColourMapType::MonochromeRed);
    qvh->set_parent (v.get_id());
    qvh->quiver_length_gain = 4.0f;
    qvh->quiver_thickness_gain = 0.005f;
    qvh->addLabel (std::format ("Vert FOV: {:.2f}{}",
                                (oces_reader.vert_fov * sm::mathconst<float>::rad2deg),
                                mplot::unicode::toUtf8(mplot::unicode::degreesign)),
                   sm::vec<float>{-4,3,0}, mplot::TextFeatures(0.12f));
    qvh->finalize();
    auto qvhp_v = v.addVisualModel (qvh);
    qvhp_v->setHide (v.view_options.test (viewopts::show_proj_fov) == false);

    auto ep = make_eye_model (v, oces_reader, ommatidia.get(), &ommatidiaColours,
                              projstr, psrad, pscentre, psr, psrax, twod_shift, nullptr);

    auto aar = oces_reader.acceptance_angle.range();
    std::cout << "Acceptance angle range: "
              << (aar.min * sm::mathconst<float>::rad2deg) << " - " << (aar.max * sm::mathconst<float>::rad2deg)
              << ", mean: " << (oces_reader.acceptance_angle.mean() * sm::mathconst<float>::rad2deg) << std::endl;
    std::cout << "Lens diameter range: " << oces_reader.diameter.range()
              << ", mean: " << oces_reader.diameter.mean() << std::endl;

    sm::flags<viewopts> last_view_options = v.view_options;
    while (!v.readyToFinish()) {
        v.waitevents (0.016);
        if (v.view_options.get() != last_view_options.get()) {
            // Projected quivers and the colourbar are independent VisualModels, so setHide():
            qvhp_h->setHide (v.view_options.test (viewopts::show_proj_fov) == false);
            qvhp_v->setHide (v.view_options.test (viewopts::show_proj_fov) == false);
            cbvp->setHide (v.view_options.test (viewopts::show_fov) == true);
            oces_origin->setHide (v.view_options.test (viewopts::show_origin) == false);
            // Everything else is part of eyevm.
            ep = make_eye_model (v, oces_reader, ommatidia.get(), &ommatidiaColours,
                                 projstr, psrad, pscentre, psr, psrax, twod_shift, ep);
            last_view_options = v.view_options;
        }
        v.render();
    }

    return 0;
}
