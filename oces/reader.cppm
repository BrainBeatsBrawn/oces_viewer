// -*- C++ -*-
module;

#include <iostream>
#include <cstdint>
#include <vector>
#include <unordered_map>
#include <stdexcept>
#include <string>
#include <sstream>

#define TINYGLTF_IMPLEMENTATION
#define STB_IMAGE_IMPLEMENTATION
#define STB_IMAGE_WRITE_IMPLEMENTATION
#if defined( WIN32 )
#pragma warning( push )
#pragma warning( disable : 4267 )
#endif
#include <tinygltf/tiny_gltf.h>
#if defined( WIN32 )
#pragma warning( pop )
#endif

export module oces.reader;

export import sm.mathconst;
export import sm.vvec;
export import sm.vec;
export import sm.mat;
export import sm.hexgrid;
import sm.hexgrid.hdf;
export import sm.bezcurvepath;
import sm.bezcurve;
import sm.centroid;
import sm.algo;

export import oces.eye;

export namespace oces
{
    // An oces::reader may compute a hex equivalent version of the eye read from file. This struct
    // holds the logic for converting from an eye of arbitrarily placed ommatidia to a hexagonally
    // arranged eye.
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
        hexeye (const oces::eye& refeye, const float iod_mult = 1.0f) { this->init (refeye, iod_mult); }

        // Compute the average height, orientation, and other parameters from any point within one
        // hex-spacing (sm::hexgrid::d) of the current hex
        void resample_eye_by_averaging (const sm::vvec<F>& eye_z, sm::vvec<F>& hex_z,
                                        const oces::eye& ref_eye, oces::eye& hex_eye,
                                        const sm::vvec<sm::vec<F, 2>>& coords2)
        {
            // Resize output containers
            hex_z.resize (this->hg.num(), F{0});
            hex_eye.orientation.resize (this->hg.num(), sm::vec<float>{});
            hex_eye.focal_offset.resize (this->hg.num(), 0.0f);
            hex_eye.diameter.resize (this->hg.num(), 0.0f);
            hex_eye.acceptance_angle.resize (this->hg.num(), 0.0f);

            const F d_thresh = this->hg.d;

            // For each hex, find all the close real data and copy the z values/average.
            for (std::uint32_t xi = 0u; xi < this->hg.num(); ++xi) {

                const sm::vec<F, 2> hp = { this->hg.d_x[xi], this->hg.d_y[xi] };

                sm::vvec<std::uint32_t> nearby = {};

                for (std::uint32_t i = 0u; i < coords2.size(); ++i) {
                    const F _d = (hp - coords2[i]).length(); // distance between coord and our hex
                    if (_d < d_thresh) { nearby.push_back (i); }
                }

                // Found nearby points, now use them
                F z_sum = F{0};
                sm::vec<float> o_sum = {};
                float f_sum = 0.0f;
                float d_sum = 0.0f;
                float a_sum = 0.0f;
                for (std::uint32_t i = 0u; i < nearby.size(); ++i) {
                    z_sum += eye_z[nearby[i]];
                    o_sum += ref_eye.orientation[nearby[i]];
                    f_sum += ref_eye.focal_offset[nearby[i]];
                    d_sum += ref_eye.diameter[nearby[i]];
                    a_sum += ref_eye.acceptance_angle[nearby[i]];
                }
                hex_z[xi] = z_sum / nearby.size();
                hex_eye.orientation[xi] = o_sum / nearby.size();
                hex_eye.focal_offset[xi] = f_sum / nearby.size();
                hex_eye.diameter[xi] = d_sum / nearby.size();
                hex_eye.acceptance_angle[xi] = a_sum / nearby.size();
            }
        }

        // @refeye: the reference OCES eye from which we are created.  @iod_mult: How many multiples
        // of interommatidial distance should we use to make the hex grid? 1 gives right
        // interommatidial spacing, but usually results in fewer ommatidia than the real eye
        void init (const oces::eye& refeye, const float iod_mult = 1.0f)
        {
            if (refeye.ready == false) {
                std::cerr << "Can't use that reference eye as it is not ready. Returning.\n";
                return;
            }

            // Copy to this eye first
            this->eye = refeye;
            this->eye.ready = false;

            // Set up hex grid. Need mean ommatidial neighbour distance to create the hexgrid, along with the approximate span.
            this->hg.init (refeye.d_mean * iod_mult, refeye.d_max * 2.0f);

            // Now use offset and angle to make a transform for the hexgrid
            F ang = refeye.central_neighbour_angles.min() + sm::mathconst<float>::pi_over_6;
            sm::mat<F, 4> tfm (sm::quaternion<F>(sm::vec<F>::uz(), ang));
            sm::vec<F> c = refeye.eye_plane_coordinates[refeye.central_omm]; // need eye_plane_orientations
            c[2] = 0;
            tfm.pretranslate (c);
            this->hg.transform (tfm);

            // Need 2D points for graham scan to find the convex hull of eye_plane_coordinates
            sm::vvec<sm::vec<F, 2>> coords2 (refeye.get_omm_per_eye());
            for (std::uint32_t i = 0; i < refeye.get_omm_per_eye(); ++i) {
                coords2[i] = refeye.eye_plane_coordinates[i].less_one_dim().template as<F>();
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

            this->resample_eye_by_averaging (eye_z, this->hg_z, refeye, this->eye, coords2);

            // Build the oces::eye, which means transforming this->hg_z, etc
            this->eye.position.resize (this->hg.num());
            auto tfm_position = refeye.eye_plane.eye_to_oces;
            for (std::uint32_t i = 0; i < this->hg.num(); ++i) {
                const sm::vec<float> p = {
                    static_cast<float>(this->hg.d_x[i]),
                    static_cast<float>(this->hg.d_y[i]),
                    static_cast<float>(this->hg_z[i])
                };
                this->eye.position[i] = (tfm_position * p).less_one_dim();
            }

            // Having built one eye, need to generate the mirror.
            this->eye.construct_mirror_eye();

            this->eye.postprocess();
        }

        // Save the eye and the hexgrid to files.
        void save (const std::string& filename_base)
        {
            // I've been using the suffix .eye in craysim programs. The format is a comma separated table.
            std::string fname_csv = filename_base + ".heye";
            std::string fname_hdf = filename_base + ".h5";
            this->eye.output_compound_ray_csv (fname_csv);
            std::cout << "Saved hexeye::eye to " << fname_csv << "\n";
            sm::hexgrid_save (this->hg, fname_hdf);
            std::cout << "Saved hexeye::hg to " << fname_hdf << "\n";
        }
    };

    struct reader
    {
        std::string filename;
        std::string base_dir = "";

        // Original eye as read from the hexgrid
        oces::eye eye;

        // Hex-equivalent eye and its associated hexgrid
        oces::hexeye<float> heye;

        // Set true if this->eye has been replaced with the hex equivalent eye (via setup_hexeye())
        bool eye_is_hex = false;

        // Get a pointer to the main eye
        oces::eye* get_eye()
        {
            if (this->eye_is_hex == true) {
                return &this->heye.eye;
            } else {
                return &this->eye;
            }
        }

        // set true to ignore any mirrors while reading
        bool ignore_mirrors = false;

        // Set true after the OCES file has been read
        bool read_success = false;

        reader() {}

        reader (const std::string& _filename, const bool& _ignore_mirrors = false)
        {
            this->filename = _filename;
            this->ignore_mirrors = _ignore_mirrors;
            this->read();
            this->postprocess();
        }

        void read (const std::string& _filename)
        {
            this->filename = _filename;
            this->read();
            this->postprocess();
        }

        // Stats etc to be computed after reading
        void postprocess() { this->eye.postprocess(); }

        // After reading the eye, you can call this to set up the hexeye
        void setup_hexeye (const float iod_mult = 1.0f)
        {
            if (this->read_success == false) {
                std::cerr << "OCES eye has not been read, so can't setup hex-equivalent. Returning.\n";
                return;
            }

            // Optional hex-equivalent eye
            heye.init (this->eye, iod_mult);
            if (heye.eye.ready == false) {
                std::cerr << "OCES eye did not convert to a hex-equivalent. Returning.\n";
                return;
            }

            this->eye_is_hex = true;
        }

        void read()
        {
            tinygltf::Model model;
            tinygltf::TinyGLTF loader;
            std::string err = "";
            std::string warn = "";

            bool r_loader = loader.LoadASCIIFromFile (&model, &err, &warn, this->filename);
            if (!warn.empty()) { std::cerr << "glTF WARNING: " << warn << std::endl; }
            if (r_loader == false) {
                std::cerr << "Failed to load GLTF file '" << filename << std::endl;
                return;
            }

            // We want to access extensions for the OCES stuff
            loader.SetStoreOriginalJSONForExtrasAndExtensions (true);

            // Calculate and store the path to the file bar the file iteself for relative includes
            this->base_dir = "";
            std::size_t slashPos = filename.find_last_of ("/\\") + 1; // (+1 to include the slash)
            if (slashPos != std::string::npos) { this->base_dir = filename.substr (0, slashPos); }

            bool oces_eyes_used = false;
            for (const auto& eu : model.extensionsUsed) {
                if (eu == "OCES_eyes") { oces_eyes_used = true; }
            }
            if (!oces_eyes_used) {
                std::cerr << "Warning: Did not find \"OCES_eyes\" in \"extensions\" section of glTF. Carrying on anyway...\n";
            }

            std::vector<std::unordered_map<std::string, int>> eyes_OmmatidialProperties;
            sm::vvec<int> ommatidialAccessors;

            // model.extensions is a map<string, Value>
            for (const auto& e : model.extensions) {

                // e.first is string, e.second is a tinygltf Value
                if (e.first != "OCES_eyes") { continue; } // We only have eyes for one extension
                auto ev = e.second;

                if (!ev.IsObject()) { continue; }
                if (!ev.Has ("ommatidialProperties")) { continue; }
                if (!ev.Has ("eyes")) { continue; }

                auto ommatidialProperties = ev.Get ("ommatidialProperties");
                if (!ommatidialProperties.IsArray()) {
                    std::cerr << "This ommatidialProperties is not an array; try next extension\n";
                    continue;
                }

                // Load ommatidialProperties
                ommatidialAccessors.resize (ommatidialProperties.Size(), 0);
                for (size_t i = 0; i < ommatidialProperties.Size(); ++i) {
                    auto ommprop = ommatidialProperties.Get(i);

                    if (ommprop.IsObject()) {
                        if (ommprop.Has("type")) {
                            auto ot = ommprop.Get("type");
                            if (ot.IsString() && ot.Get<std::string>() == "ACCESSOR") {
                                if (ommprop.Has("value")) {
                                    auto ov = ommprop.Get("value");
                                    if (ov.IsInt()) {
                                        // ommatidialProperty[i] is ACCESSOR with value ov.Get<int>()
                                        ommatidialAccessors[i] = ov.Get<int>();
                                    }
                                }
                            }
                        }
                    }
                }

                // Read "eyes" from "OCES_eyes"
                auto eyes = ev.Get ("eyes");
                if (eyes.IsArray()) {

                    eyes_OmmatidialProperties.resize (eyes.Size());

                    for (size_t i = 0; i < eyes.Size(); ++i) {
                        // Process eye
                        if (eyes.Get(i).Get("type").Get<std::string>() != "POINT_OMMATIDIAL") {
                            std::cerr << "Don't know how to process an OCES eye of type '"
                                      << eyes.Get(i).Get("type").Get<std::string>() << "'\n";
                            continue;
                        }

                        auto op = eyes.Get(i).Get("ommatidialProperties");

                        if (!op.IsObject()) {
                            std::cerr << "Badly formed OCES glTF (OCES_eyes.ommatidialProperties is not a JSON object)\n";
                            continue;
                        }
                        if (!(op.Has("POSITION") && op.Has("ORIENTATION") && op.Has("FOCAL_OFFSET") && op.Has("DIAMETER"))) {
                            std::cerr << "Badly formed OCES glTF (OCES_eyes.ommatidialProperties is not a JSON object)\n";
                            continue;
                        }

                        std::cerr << "Processing eye " << eyes.Get(i).Get("name").Get<std::string>()
                                  << " of type " << eyes.Get(i).Get("type").Get<std::string>() << std::endl;

                        // Good to go
                        eyes_OmmatidialProperties[i]["POSITION"] = op.Get("POSITION").Get<int>();
                        eyes_OmmatidialProperties[i]["ORIENTATION"] = op.Get("ORIENTATION").Get<int>();
                        eyes_OmmatidialProperties[i]["FOCAL_OFFSET"] = op.Get("FOCAL_OFFSET").Get<int>();
                        eyes_OmmatidialProperties[i]["DIAMETER"] = op.Get("DIAMETER").Get<int>();
                    }
                }

                // Can now read the buffers into our member attributes position, orientation, etc
                for (auto eye : eyes_OmmatidialProperties) {
                    int sz = static_cast<int>(ommatidialAccessors.size());
                    if (eye["POSITION"] < sz) { this->get_buffer (model, ommatidialAccessors[eye["POSITION"]], this->eye.position); }
                    if (eye["ORIENTATION"] < sz) { this->get_buffer (model, ommatidialAccessors[eye["ORIENTATION"]], this->eye.orientation); }
                    if (eye["FOCAL_OFFSET"] < sz) { this->get_buffer (model, ommatidialAccessors[eye["FOCAL_OFFSET"]], this->eye.focal_offset); }
                    if (eye["DIAMETER"] < sz) { this->get_buffer (model, ommatidialAccessors[eye["DIAMETER"]], this->eye.diameter); }
                }

                // Compound-ray eye files use acceptance angle, rather than optical lens diameter
                this->eye.acceptance_angle.resize (this->eye.diameter.size(), 0.0f);
                for (size_t i = 0; i < this->eye.diameter.size(); i++) {
                    this->eye.acceptance_angle[i] = 2.0f * std::atan2 (this->eye.diameter[i] / 2.0f, std::abs (this->eye.focal_offset[i]));
                }

                // Read "mirrorPlanes" from "OCES_eyes"
                if (ev.Has("mirrorPlanes")) {
                    auto mplanes = ev.Get("mirrorPlanes");
                    if (mplanes.IsArray()) {
                        for (size_t i = 0; i < mplanes.Size(); ++i) {
                            oces::plane mp;
                            if (mplanes.Get(i).Get("normal").Get<std::string>() == "X"
                                || mplanes.Get(i).Get("normal").Get<std::string>() == "FRONTAL"
                                || mplanes.Get(i).Get("normal").Get<std::string>() == "LEFT"
                                || mplanes.Get(i).Get("normal").Get<std::string>() == "RIGHT") {
                                mp.normal = sm::vec<float>::ux();
                            } else if (mplanes.Get(i).Get("normal").Get<std::string>() == "Y"
                                       || mplanes.Get(i).Get("normal").Get<std::string>() == "TRANSVERSE"
                                       || mplanes.Get(i).Get("normal").Get<std::string>() == "DORSOVENTRAL"
                                       || mplanes.Get(i).Get("normal").Get<std::string>() == "UP"
                                       || mplanes.Get(i).Get("normal").Get<std::string>() == "DOWN") {
                                mp.normal = sm::vec<float>::uy();
                            } else if (mplanes.Get(i).Get("normal").Get<std::string>() == "Z"
                                       || mplanes.Get(i).Get("normal").Get<std::string>() == "SAGITTAL"
                                       || mplanes.Get(i).Get("normal").Get<std::string>() == "ANTERIOR"
                                       || mplanes.Get(i).Get("normal").Get<std::string>() == "POSTERIOR"
                                       || mplanes.Get(i).Get("normal").Get<std::string>() == "ANTEPOSTERIOR"
                                       || mplanes.Get(i).Get("normal").Get<std::string>() == "FORWARD"
                                       || mplanes.Get(i).Get("normal").Get<std::string>() == "BACKWARD"
                                       || mplanes.Get(i).Get("normal").Get<std::string>() == "FRONT"
                                       || mplanes.Get(i).Get("normal").Get<std::string>() == "BACK") {
                                mp.normal = sm::vec<float>::uz();
                            } else {
                                // could also be an array of 3 numbers to specify a normal vector
                                if (mplanes.Get(i).Get("normal").IsArray()) {
                                    auto nrm = mplanes.Get(i).Get("normal");
                                    if (nrm.Size() == 3) {
                                        for (size_t j = 0; j < 3; ++j) {
                                            double element = nrm.Get(j).Get<double>();
                                            mp.normal[j] = static_cast<float>(element);
                                        }
                                    } else {
                                        throw std::runtime_error ("mirror plane normal specified with wrong number of dimensions");
                                    }
                                }
                            }

                            // may also need code to get the mirrorplane position
                            if (mplanes.Get(i).Get("position").IsArray()) {
                                auto posn = mplanes.Get(i).Get("position");
                                if (posn.Size() == 3) {
                                    for (size_t j = 0; j < 3; ++j) {
                                        double element = posn.Get(j).Get<double>();
                                        mp.origin[j] = static_cast<float>(element);
                                    }
                                } else {
                                    throw std::runtime_error ("mirror plane position specified with wrong number of dimensions");
                                }
                            }

                            this->eye.mirrorplanes.push_back (mp);
                        }
                    }
                }

                // Compute statistics on the eye
                std::cerr << "Number of ommatidia: " << this->eye.position.size() << std::endl;
                std::cerr << "Optical diameter mean/std: "
                          << this->eye.diameter.mean() << " (" << this->eye.diameter.std() << ")\n";
                std::cerr << "Acceptance angle mean/std degrees: "
                          << this->eye.acceptance_angle.mean() * sm::mathconst<float>::rad2deg
                          << " (" << this->eye.acceptance_angle.std() * sm::mathconst<float>::rad2deg << ")\n";
                // FOV
                // Find mean direction
                sm::vec<float> mean_dir = {};
                for (auto ori : this->eye.orientation) {
                    sm::vec<float> orient = ori;
                    orient.renormalize();
                    mean_dir += orient;
                }
                mean_dir /= this->eye.orientation.size();
                sm::vvec<float> ang_from_mean;
                for (auto ori : this->eye.orientation) {
                    auto a = ori.angle (mean_dir);
                    ang_from_mean.push_back (a);
                }
                std::cerr << "ang_from_mean max " << ang_from_mean.max() * sm::mathconst<float>::rad2deg << std::endl;

                // Act on mirror planes and add to position arrays
                if (!this->eye.mirrorplanes.empty() && !this->ignore_mirrors) {
                    this->eye.construct_mirror_eye();
                }
            }

            // Process nodes (to get the head mesh)
            sm::mat<float, 4> root_transform;
            std::vector<int32_t> root_nodes (model.nodes.size(), 1);
            for (auto& gltf_node : model.nodes) {
                for (int32_t child : gltf_node.children) { root_nodes[child] = 0; }
            }

            for (size_t i = 0; i < root_nodes.size(); ++i) {
                if (!root_nodes[i]) { continue; }
                auto& gltf_node = model.nodes[i];
                this->process_node (model, gltf_node, root_transform);
            }

            this->read_success = true;
        }

        /**
         * Copy data from the buffer identified by accessor_idx into the output vvec.
         *
         * \tparam T The type of the data elements in output. May be a scalar such as float or double, or a
         * fixed size data type such as sm::vec<float, 3> or float3
         *
         * \param model The initialized (by a tinygltf::TinyGLTF loader) tinygltf model reference.
         *
         * \param accessor_idx An integer index for the accessor to the glTF buffer
         *
         * \param output The output vvec. This will be resized, then filled with a *copy* of the data held
         * in the TinyGLTF model.
         */
        template<typename T>
        void get_buffer (const tinygltf::Model& model, const int32_t accessor_idx, std::vector<T>& output)
        {
            constexpr bool debug_get_buffer = false;

            if (accessor_idx == -1) { return; }

            const tinygltf::Accessor& gltf_accessor      = model.accessors[accessor_idx];
            const tinygltf::BufferView& gltf_buffer_view = model.bufferViews[gltf_accessor.bufferView];
            const int32_t elmt_cmpt_byte_size = tinygltf::GetComponentSizeInBytes(gltf_accessor.componentType);
            const int32_t cmpts_in_type  = tinygltf::GetNumComponentsInType(gltf_accessor.type);

            if (elmt_cmpt_byte_size == (static_cast<int32_t>(sizeof(T)) / cmpts_in_type)) {
                output.resize (gltf_accessor.count);
                // copy data from model.buffers[gltf_buffer_view.buffer].data (vector<unsigned char>) to vvec_of_vec.
                if constexpr (debug_get_buffer) {
                    std::cerr << "Memcpy " << gltf_accessor.count * elmt_cmpt_byte_size * cmpts_in_type
                              << " bytes from accessor index " << accessor_idx << ", buffer view byte offset is " << gltf_buffer_view.byteOffset << "\n";
                }
                std::memcpy (output.data(),
                             model.buffers[gltf_buffer_view.buffer].data.data() + gltf_buffer_view.byteOffset,
                             gltf_accessor.count * elmt_cmpt_byte_size * cmpts_in_type);
            } else if (elmt_cmpt_byte_size < (static_cast<int32_t>(sizeof(T)) / cmpts_in_type)) {
                std::stringstream ee;
                ee << "Failed to memcpy in get_buffer because elmt_cmpt_byte_size " << elmt_cmpt_byte_size
                   << " < sizeof(T) = " << static_cast<int32_t>(sizeof(T))
                   << "\n sizeof(T)/cmpts_in_type = " << (static_cast<int32_t>(sizeof(T)) / cmpts_in_type);
                throw std::runtime_error (ee.str());
            } else {
                std::stringstream ee;
                ee << "Failed to memcpy in get_buffer!\n cmpts_in_type = " << cmpts_in_type
                   << "\n elmt_cmpt_byte_size = " <<  elmt_cmpt_byte_size
                   << "\n sizeof(T) = " << static_cast<int32_t>(sizeof(T))
                   << "\n sizeof(T)/cmpts_in_type = " << (static_cast<int32_t>(sizeof(T)) / cmpts_in_type);
                throw std::runtime_error (ee.str());
            }
        }

        // Process nodes to read head mesh
        void process_node (const tinygltf::Model& model, const tinygltf::Node& gltf_node, const sm::mat<float, 4>& parent_matrix)
        {
            constexpr bool debug_gltf = true;

            std::cerr << "Process node " << gltf_node.name << std::endl;

            sm::mat<float, 4> translation;
            if (!gltf_node.translation.empty()) {
                auto tr = sm::vec<double>{ gltf_node.translation[0], gltf_node.translation[1], gltf_node.translation[2] };
                sm::vec<float> trf = tr.as<float>();
                translation.translate (trf);
            }
            sm::mat<float, 4> rotation;
            if (!gltf_node.rotation.empty()) {
                sm::quaternion<double> q(gltf_node.rotation[3], gltf_node.rotation[0], gltf_node.rotation[1], gltf_node.rotation[2]);
                rotation.rotate (q);
            }
            sm::mat<float, 4> scale;
            if (!gltf_node.scale.empty()) {
                scale.scale (sm::vec<double>({ gltf_node.scale[0], gltf_node.scale[1], gltf_node.scale[2] }).as<float>());
            }
            sm::mat<float, 4> matrix;
            if (!gltf_node.matrix.empty()) {
                for (uint32_t i = 0; i < 16; ++i) { matrix.arr[i] = static_cast<float>(gltf_node.matrix[i]); }
            }
            const sm::mat<float, 4> node_xform = parent_matrix * matrix * translation * rotation * scale;

            if (gltf_node.mesh != -1) {
                const auto& gltf_mesh = model.meshes[gltf_node.mesh];
                if constexpr (debug_gltf == true) {
                    std::cerr << "Processing glTF mesh: '" << gltf_mesh.name << "'\n";
                    std::cerr << "\tNum mesh primitive groups: " << gltf_mesh.primitives.size() << std::endl;
                }
                if (gltf_node.name == "head") {
                    std::cerr << "Have head node; read mesh(es)\n";
                    for (auto& gltf_primitive : gltf_mesh.primitives) {
                        if (gltf_primitive.mode != TINYGLTF_MODE_TRIANGLES) { throw std::runtime_error ("Non-triangle primitive"); }
                        std::cerr << "Have head mesh; reading it\n";
                        this->eye.head_mesh.name = gltf_mesh.name;
                        // Indices
                        {
                            const tinygltf::Accessor& _accessor = model.accessors[gltf_primitive.indices];
                            const int32_t elmt_cmpt_byte_size = tinygltf::GetComponentSizeInBytes (_accessor.componentType);
                            const int32_t cmpts_in_type  = tinygltf::GetNumComponentsInType (_accessor.type);
                            if (cmpts_in_type != 1) { throw std::runtime_error ("Expect 1 component in type for indices"); }
                            if (elmt_cmpt_byte_size == 4) {
                                this->get_buffer<uint32_t> (model, gltf_primitive.indices, this->eye.head_mesh.indices);
                            } else if (elmt_cmpt_byte_size == 2) {
                                std::vector<uint16_t> hmi16;
                                this->get_buffer<uint16_t> (model, gltf_primitive.indices, hmi16);
                                // copy to head_mesh.indices
                                this->eye.head_mesh.indices.resize (hmi16.size());
                                for (uint32_t i = 0; i < hmi16.size(); ++i) { this->eye.head_mesh.indices[i] = hmi16[i]; }
                            } else {
                                throw std::runtime_error ("Deal with 1 byte or 8 byte index size");
                            }
                        }
                        this->eye.head_mesh.transform = node_xform;
                        if constexpr (debug_gltf == true) { std::cerr << "\t\tNum triangles is indices size/3: " << this->eye.head_mesh.indices.size() / 3 << std::endl; }
                        // Positions
                        assert (gltf_primitive.attributes.find( "POSITION" ) !=  gltf_primitive.attributes.end());
                        const int32_t pos_accessor_idx = gltf_primitive.attributes.at ("POSITION");
                        this->get_buffer<sm::vec<float>> (model, pos_accessor_idx, this->eye.head_mesh.positions);
                        if constexpr (debug_gltf == true) { std::cerr << "\t\tNum vertices(positions count): " << this->eye.head_mesh.positions.size() << std::endl; }

                        const auto& pos_gltf_accessor = model.accessors[pos_accessor_idx];
                        sm::vec<double> minvals = { pos_gltf_accessor.minValues[0], pos_gltf_accessor.minValues[1], pos_gltf_accessor.minValues[2] };
                        sm::vec<double> maxvals = { pos_gltf_accessor.maxValues[0], pos_gltf_accessor.maxValues[1], pos_gltf_accessor.maxValues[2] };
                        this->eye.head_mesh.object_aabb = sm::interval<sm::vec<float>> (minvals.as<float>(), maxvals.as<float>());
                        this->eye.head_mesh.world_aabb = sm::interval<sm::vec<float>> ((node_xform * this->eye.head_mesh.object_aabb.min).less_one_dim(),
                                                                                       (node_xform * this->eye.head_mesh.object_aabb.max).less_one_dim());
                        // Normals
                        auto normal_accessor_iter = gltf_primitive.attributes.find ("NORMAL");
                        if (normal_accessor_iter != gltf_primitive.attributes.end()) {
                            if constexpr (debug_gltf == true) { std::cerr << "\t\tHas vertex normals: true\n"; }
                            const int32_t normal_accessor_idx = gltf_primitive.attributes.at ("NORMAL");
                            this->get_buffer<sm::vec<float>> (model, normal_accessor_idx, this->eye.head_mesh.normals);
                        }
                        // omit to read texture coordinates
                        // omit handling of vertex colours
                    }
                }

            } else if (!gltf_node.children.empty()) {
                for (int32_t child : gltf_node.children) {
                    this->process_node (model, model.nodes[child], node_xform);
                }
            }
        }
    };

} // namespace oces
