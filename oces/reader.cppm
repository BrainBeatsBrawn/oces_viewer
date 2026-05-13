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
export import sm.geometry;

export namespace oces
{
    struct mirrorplane
    {
        sm::vec<float> position = {};
        sm::vec<float> normal = {1, 0, 0};
    };

    // This struct matches the layout of mplot::meshgroup
    struct meshgroup
    {
        std::string name;
        sm::mat<float, 4> transform;
        sm::vvec<uint32_t> indices;
        sm::vvec<sm::vec<float>> positions;
        sm::vvec<sm::vec<float>> normals;
        sm::vvec<sm::vec<float>> colours;
        sm::range<sm::vec<float>> object_aabb;
        sm::range<sm::vec<float>> world_aabb;
        std::array<float, 3> single_colour = {0};
    };

    struct reader
    {
        std::string filename;
        std::string base_dir = "";

        // These probably want to go into an 'eye' struct
        sm::vvec<sm::vec<float, 3>> position = {};    // Units: m
        sm::vvec<sm::vec<float, 3>> orientation = {}; // Units: m
        sm::vvec<float> focal_offset = {};            // Units: m
        sm::vvec<float> diameter = {};                // Optical diameter. Units: m
        sm::vvec<float> acceptance_angle = {};        // Derived from diameter and focal offset. Units: radians.

        // Horizontal field of view (about the up axis, which is y in OCES)
        float horz_fov = 0.0f;
        // Vertical field of view (about the fwds axis, which is z in OCES)
        float vert_fov = 0.0f;

        // Horizontally projected orientations, for visualization
        sm::vvec<sm::vec<float, 3>> h_plane_orientation;
        sm::vvec<sm::vec<float, 3>> h_plane_position;
        // Vertically projected orientations, for visualization
        sm::vvec<sm::vec<float, 3>> v_plane_orientation;
        sm::vvec<sm::vec<float, 3>> v_plane_position;

        // A second eye may be mirrored by a mirrorplane
        std::vector<oces::mirrorplane> mirrorplanes;
        // Store the matrix for the mirrors defined in mirrorplanes
        std::vector<sm::mat<float, 4>> mirrors;
        // set true to ignore any mirrors
        bool ignore_mirrors = false;

        // If we have a head mesh, store it here
        oces::meshgroup head_mesh;

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
        }

        void compute_fov_max()
        {
            auto horz_fov_r = sm::range<float>::search_initialized();
            auto vert_fov_r = sm::range<float>::search_initialized();

            std::uint32_t omm_per_eye = this->orientation.size();
            if (!this->mirrorplanes.empty() && !this->ignore_mirrors) {
                // Two eyes are stored in this->orientation.
                omm_per_eye /= 2u;
            }

            this->h_plane_orientation.resize (omm_per_eye);
            this->v_plane_orientation.resize (omm_per_eye);
            this->h_plane_position.resize (omm_per_eye);
            this->v_plane_position.resize (omm_per_eye);

            for (std::uint32_t i = 0; i < omm_per_eye; ++i) {
                // horz fov
                auto onto_y_i = sm::geometry::vector_plane_projection (sm::vec<>::uy(), this->orientation[i]);
                this->h_plane_orientation[i] = onto_y_i;
                this->h_plane_position[i] = sm::geometry::vector_plane_projection (sm::vec<>::uy(), this->position[i]);
                // vert fov
                auto onto_z_i = sm::geometry::vector_plane_projection (sm::vec<>::uz(), this->orientation[i]);
                this->v_plane_orientation[i] = onto_z_i;
                this->v_plane_position[i] = sm::geometry::vector_plane_projection (sm::vec<>::uz(), this->position[i]);

                // Compare angle to all others
                for (std::uint32_t j = 0; j < omm_per_eye; ++j) {
                    if (j != i) {
                        // project onto the plane
                        auto onto_y_j = sm::geometry::vector_plane_projection (sm::vec<>::uy(), this->orientation[j]); // horz fov
                        auto onto_z_j = sm::geometry::vector_plane_projection (sm::vec<>::uz(), this->orientation[j]); // vert fov

                        float horz_ang = onto_y_i.angle (onto_y_j);
                        float vert_ang = onto_z_i.angle (onto_z_j);
                        horz_fov_r.update (horz_ang);
                        vert_fov_r.update (vert_ang);
                    }
                }
            }

            this->horz_fov = horz_fov_r.max;
            this->vert_fov = vert_fov_r.max;
        }

        // Stats etc to be computed after reading
        void postprocess()
        {
            this->compute_fov_max();
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
                std::cout << "Warning: Did not find \"OCES_eyes\" in \"extensions\" section of glTF. Carrying on anyway...\n";
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
                    std::cout << "This ommatidialProperties is not an array; try next extension\n";
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
                            std::cout << "Don't know how to process an OCES eye of type '"
                                      << eyes.Get(i).Get("type").Get<std::string>() << "'\n";
                            continue;
                        }

                        auto op = eyes.Get(i).Get("ommatidialProperties");

                        if (!op.IsObject()) {
                            std::cout << "Badly formed OCES glTF (OCES_eyes.ommatidialProperties is not a JSON object)\n";
                            continue;
                        }
                        if (!(op.Has("POSITION") && op.Has("ORIENTATION") && op.Has("FOCAL_OFFSET") && op.Has("DIAMETER"))) {
                            std::cout << "Badly formed OCES glTF (OCES_eyes.ommatidialProperties is not a JSON object)\n";
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
                    if (eye["POSITION"] < sz) { this->get_buffer (model, ommatidialAccessors[eye["POSITION"]], this->position); }
                    if (eye["ORIENTATION"] < sz) { this->get_buffer (model, ommatidialAccessors[eye["ORIENTATION"]], this->orientation); }
                    if (eye["FOCAL_OFFSET"] < sz) { this->get_buffer (model, ommatidialAccessors[eye["FOCAL_OFFSET"]], this->focal_offset); }
                    if (eye["DIAMETER"] < sz) { this->get_buffer (model, ommatidialAccessors[eye["DIAMETER"]], this->diameter); }
                }

                // Compound-ray eye files use acceptance angle, rather than optical lens diameter
                this->acceptance_angle.resize (this->diameter.size(), 0.0f);
                for (size_t i = 0; i < this->diameter.size(); i++) {
                    this->acceptance_angle[i] = 2.0f * std::atan2 (this->diameter[i] / 2.0f, std::abs (this->focal_offset[i]));
                }

                // Read "mirrorPlanes" from "OCES_eyes"
                if (ev.Has("mirrorPlanes")) {
                    auto mplanes = ev.Get("mirrorPlanes");
                    if (mplanes.IsArray()) {
                        for (size_t i = 0; i < mplanes.Size(); ++i) {
                            oces::mirrorplane mp;
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
                                        mp.position[j] = static_cast<float>(element);
                                    }
                                } else {
                                    throw std::runtime_error ("mirror plane position specified with wrong number of dimensions");
                                }
                            }

                            this->mirrorplanes.push_back (mp);
                        }
                    }
                }

                // Compute statistics on the eye
                std::cerr << "Number of ommatidia: " << this->position.size() << std::endl;
                std::cerr << "Optical diameter mean/std: "
                          << this->diameter.mean() << " (" << this->diameter.std() << ")\n";
                std::cerr << "Acceptance angle mean/std degrees: "
                          << this->acceptance_angle.mean() * sm::mathconst<float>::rad2deg
                          << " (" << this->acceptance_angle.std() * sm::mathconst<float>::rad2deg << ")\n";
                // FOV
                // Find mean direction
                sm::vec<float> mean_dir = {};
                for (auto ori : this->orientation) {
                    sm::vec<float> orient = ori;
                    orient.renormalize();
                    mean_dir += orient;
                }
                mean_dir /= this->orientation.size();
                sm::vvec<float> ang_from_mean;
                for (auto ori : this->orientation) {
                    auto a = ori.angle (mean_dir);
                    ang_from_mean.push_back (a);
                }
                std::cerr << "ang_from_mean max " << ang_from_mean.max() * sm::mathconst<float>::rad2deg << std::endl;

                // Act on mirror planes and add to position arrays
                if (!this->mirrorplanes.empty() && !this->ignore_mirrors) {

                    // Get size for one eye
                    size_t sz = this->position.size();

                    // Make space
                    this->position.resize (2 * sz);
                    this->orientation.resize (2 * sz);
                    this->focal_offset.resize (2 * sz);
                    this->diameter.resize (2 * sz);
                    this->acceptance_angle.resize (2 * sz);

                    sm::mat<float, 4> mirror = sm::mat<float, 4>::reflection (this->mirrorplanes[0].position, this->mirrorplanes[0].normal);
                    this->mirrors.push_back (mirror); // Saved for client code to use
                    for (size_t i = 0; i < sz; ++i) {
                        // Mirror position and direction
                        sm::vec<float> mpos = (mirror * this->position[i]).less_one_dim();
                        sm::vec<float> mdir = (mirror * this->orientation[i]).less_one_dim();
                        position[sz + i] = mpos;
                        orientation[sz + i] = mdir;
                        // Focal offset, diameter and acceptance angle are simply copied
                        focal_offset[sz + i] = focal_offset[i];
                        diameter[sz + i] = diameter[i];
                        acceptance_angle[sz + i] = acceptance_angle[i];
                    }
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
                        head_mesh.name = gltf_mesh.name;
                        // Indices
                        {
                            const tinygltf::Accessor& _accessor = model.accessors[gltf_primitive.indices];
                            const int32_t elmt_cmpt_byte_size = tinygltf::GetComponentSizeInBytes (_accessor.componentType);
                            const int32_t cmpts_in_type  = tinygltf::GetNumComponentsInType (_accessor.type);
                            if (cmpts_in_type != 1) { throw std::runtime_error ("Expect 1 component in type for indices"); }
                            if (elmt_cmpt_byte_size == 4) {
                                this->get_buffer<uint32_t> (model, gltf_primitive.indices, head_mesh.indices);
                            } else if (elmt_cmpt_byte_size == 2) {
                                std::vector<uint16_t> hmi16;
                                this->get_buffer<uint16_t> (model, gltf_primitive.indices, hmi16);
                                // copy to head_mesh.indices
                                head_mesh.indices.resize (hmi16.size());
                                for (uint32_t i = 0; i < hmi16.size(); ++i) { head_mesh.indices[i] = hmi16[i]; }
                            } else {
                                throw std::runtime_error ("Deal with 1 byte or 8 byte index size");
                            }
                        }
                        head_mesh.transform = node_xform;
                        if constexpr (debug_gltf == true) { std::cerr << "\t\tNum triangles is indices size/3: " << head_mesh.indices.size() / 3 << std::endl; }
                        // Positions
                        assert (gltf_primitive.attributes.find( "POSITION" ) !=  gltf_primitive.attributes.end());
                        const int32_t pos_accessor_idx = gltf_primitive.attributes.at ("POSITION");
                        this->get_buffer<sm::vec<float>> (model, pos_accessor_idx, head_mesh.positions);
                        if constexpr (debug_gltf == true) { std::cerr << "\t\tNum vertices(positions count): " << head_mesh.positions.size() << std::endl; }

                        const auto& pos_gltf_accessor = model.accessors[pos_accessor_idx];
                        sm::vec<double> minvals = { pos_gltf_accessor.minValues[0], pos_gltf_accessor.minValues[1], pos_gltf_accessor.minValues[2] };
                        sm::vec<double> maxvals = { pos_gltf_accessor.maxValues[0], pos_gltf_accessor.maxValues[1], pos_gltf_accessor.maxValues[2] };
                        head_mesh.object_aabb = sm::range<sm::vec<float>> (minvals.as<float>(), maxvals.as<float>());
                        head_mesh.world_aabb = sm::range<sm::vec<float>> ((node_xform * head_mesh.object_aabb.min).less_one_dim(),
                                                                          (node_xform * head_mesh.object_aabb.max).less_one_dim());
                        // Normals
                        auto normal_accessor_iter = gltf_primitive.attributes.find ("NORMAL");
                        if (normal_accessor_iter != gltf_primitive.attributes.end()) {
                            if constexpr (debug_gltf == true) { std::cerr << "\t\tHas vertex normals: true\n"; }
                            const int32_t normal_accessor_idx = gltf_primitive.attributes.at ("NORMAL");
                            this->get_buffer<sm::vec<float>> (model, normal_accessor_idx, head_mesh.normals);
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

        void output_compound_ray_csv()
        {
            if (this->position.size() == this->orientation.size()
                && this->position.size() == this->focal_offset.size()
                && this->position.size() == this->acceptance_angle.size()) {
                for (size_t i = 0; i < this->position.size(); ++i) {
                    std::cout << this->position[i].str_comma_separated (' ') << " "
                              << this->orientation[i].str_comma_separated (' ')
                              << " " << this->acceptance_angle[i]
                              << " " << std::abs(this->focal_offset[i])
                              << std::endl;
                }
            } else {
                std::cerr << "position, orientation, focal_offset and diameter/acceptance_angle should all have the same number of elements.\n";
            }
        }
    };

} // namespace oces
