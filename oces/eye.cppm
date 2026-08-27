// -*- C++ -*-
module;

#include <iostream>
#include <cstdint>
#include <vector>
#include <string>
#include <fstream>

export module oces.eye;

export import sm.mathconst;
export import sm.vvec;
export import sm.vec;
export import sm.mat;
export import sm.geometry;
import sm.centroid;

export namespace oces
{
    // A 2D plane in 3D, defined by an origin and a normal
    struct plane
    {
        sm::vec<float> origin = {};
        sm::vec<float> normal = sm::vec<float>::uz();
    };

    // For an *eye* plane, we define a 'forwards' vector in the plane, which is derived from the
    // uz() direction in the OCES frame of reference. This is therefore a coordinate frame.
    struct eyeplane
    {
        oces::plane p;
        // ez is p.normal
        sm::vec<float> ex = sm::vec<float>::ux();
        sm::vec<float> ey = sm::vec<float>::uy();
        // This transform takes a coordinate in the eye frame, and converts it to a coordinate in
        // the OCES model frame.
        sm::mat<float, 4> eye_to_oces = {};
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
        sm::interval<sm::vec<float>> object_aabb;
        sm::interval<sm::vec<float>> world_aabb;
        std::array<float, 3> single_colour = {0};
    };

    // The OCES eye contains the actual data about the ommatidia in the eye (along with an optional
    // head mesh and mirrors)
    struct eye
    {
        sm::vvec<sm::vec<float, 3>> position = {};    // Units: m
        sm::vvec<sm::vec<float, 3>> orientation = {}; // Units: m
        sm::vvec<float> focal_offset = {};            // Units: m
        sm::vvec<float> diameter = {};                // Optical diameter. Units: m
        sm::vvec<float> acceptance_angle = {};        // Derived from diameter and focal offset. Units: radians.

        // position, projected onto the eye plane, and in the coordinate system defined by the eye
        // plane's normal (and what other vector? Probably a projection of uz, the forward axis in
        // the world.). Maybe just need transform * position.
        sm::vvec<sm::vec<float, 3>> eye_plane_coordinates = {};

        // The orientations, rotated into the eye plane
        sm::vvec<sm::vec<float, 3>> eye_plane_orientation = {}; // transform with this->eye_plane.eye_to_oces

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

        // Mean neighbour-neighbour distance of the closest 6 ommatidial neighbours
        float d_mean = 0.0f;

        // Most central ommatidium index
        std::uint32_t central_omm = std::numeric_limits<std::uint32_t>::max();
        // Neighbour indices
        sm::vec<sm::vec<float>, 6> central_neighbours = {};
        // The angles
        sm::vec<float, 6> central_neighbour_angles;

        // Maximum distance between any two ommatidia in the eye (just one, not a mirrored pair)
        float d_max = 0.0f;

        // A second eye may be mirrored by a mirrorplane
        std::vector<oces::plane> mirrorplanes;
        // Store the matrix for the mirrors defined in mirrorplanes
        std::vector<sm::mat<float, 4>> mirrors;

        // The first/main eye could have a plane whose normal is the average of the orientations and
        // origin is the mean location of the ommatidia. It also has an 'x' axis related to the 'z'
        // axis in the world frame. note that this plane has nothing to do with the OCES standard;
        // it's an addition by Seb to help replace the arbitrary, measured ommatidial positions with
        // an equivalent, perfectly hexagonal grid of ommatidia for simulation-friendly eyes.
        oces::eyeplane eye_plane;

        // If we have a head mesh, store it here
        oces::meshgroup head_mesh;

        // Ready to be used?
        bool ready = false;

        // Common function to find number of ommatidium in a single eye
        std::uint32_t get_omm_per_eye() const
        {
            std::uint32_t omm_per_eye = this->orientation.size();
            if (!this->mirrorplanes.empty()) {
                // Two eyes are stored in this->orientation.
                omm_per_eye /= 2u;
            }
            return omm_per_eye;
        }

        // In eye plane, find the most central OMM
        void compute_central()
        {
            const std::uint32_t omm_per_eye = this->get_omm_per_eye();

            auto ctrd = sm::algo::centroid (this->eye_plane_coordinates).less_one_dim();

            float min_d_c = std::numeric_limits<float>::max();
            for (std::uint32_t i = 0; i < omm_per_eye; ++i) {
                sm::vec<float, 2> pi = this->eye_plane_coordinates[i].less_one_dim();
                float d_c = (pi - ctrd).length();
                if (d_c < min_d_c) {
                    this->central_omm = i; // index of central ommatidium
                    min_d_c = d_c;
                }
            }

            if (central_omm > this->eye_plane_coordinates.size()) {
                std::cerr << "Out of range of position?\n";
                return;
            }

            sm::vec<float, 6> nearest_6 = {};
            nearest_6.set_from (std::numeric_limits<float>::max());
            const auto& pi = this->eye_plane_coordinates[this->central_omm];
            for (std::uint32_t j = 0; j < omm_per_eye; ++j) {
                if (j == this->central_omm) { continue; }
                const auto& pj = this->eye_plane_coordinates[j];
                const float dij = (pi - pj).length();
                // If dij is less than any, replace the biggest
                auto kbig = nearest_6.argmax();
                for (std::uint32_t k = 0; k < 6; ++k) {
                    if (dij <= nearest_6[k]) {
                        nearest_6[kbig] = dij;
                        this->central_neighbours[kbig] = pj;
                        break;
                    }
                }
            }

            for (std::uint32_t k = 0; k < 6; ++k) {
                this->central_neighbour_angles[k] = pi.less_one_dim().angle (this->central_neighbours[k].less_one_dim());
            }
        }

        // Find a mean neighbour-to-neighbour distance. (EyeVisual does something like this, too,
        // but it finds the min distance for *each* ommatidium rather than an average.
        // While at it, find d_max - the maximum straight line distance between any pair of ommatidia
        void compute_neighbour_distance()
        {
            const std::uint32_t omm_per_eye = this->get_omm_per_eye();

            if (omm_per_eye < 7) {
                // we'd not fill nearest_6 with valid distances
                std::cerr << "compute_neighbour_distance(): Too few ommatidia for this function; returning\n";
                return;
            }

            this->d_mean = 0.0f;
            this->d_max = 0.0f;
            for (std::uint32_t i = 0; i < omm_per_eye; ++i) {
                sm::vec<float, 6> nearest_6 = {};
                nearest_6.set_from (std::numeric_limits<float>::max());
                const auto& pi = this->position[i];
                for (std::uint32_t j = 0; j < omm_per_eye; ++j) {
                    if (j == i) { continue; }
                    const auto& pj = this->position[j];
                    const float dij = (pi - pj).length();
                    this->d_max = dij > this->d_max ? dij : this->d_max;
                    // If dij is less than any, replace the biggest
                    auto kbig = nearest_6.argmax();
                    for (std::uint32_t k = 0; k < 6; ++k) {
                        if (dij <= nearest_6[k]) {
                            nearest_6[kbig] = dij;
                            break;
                        }
                    }
                }
                this->d_mean += nearest_6.mean();
            }
            this->d_mean /= omm_per_eye;
        }

        // Project each ommatidium position onto the eye plane. This becomes a vector of 3D values
        // that have their first two coordinates in the plane, with the third coordinate being the
        // height above/below the plane.
        void compute_projections_on_eye_plane()
        {
            const std::uint32_t omm_per_eye = this->get_omm_per_eye();

            this->eye_plane_coordinates.resize (omm_per_eye, {});
            this->eye_plane_orientation.resize (omm_per_eye, {});

            auto tfm_position = this->eye_plane.eye_to_oces.inverse();
            auto tfm_orientation = tfm_position.linear();
            for (std::uint32_t i = 0; i < omm_per_eye; ++i) {
                this->eye_plane_coordinates[i] = (tfm_position * this->position[i]).less_one_dim();
                this->eye_plane_orientation[i] = tfm_orientation * this->orientation[i];
            }
        }

        // Find the normal and origin of this->eye_plane
        void compute_eye_plane()
        {
            const std::uint32_t omm_per_eye = this->get_omm_per_eye();
            this->eye_plane.p.origin = {};
            this->eye_plane.p.normal = {};
            for (std::uint32_t i = 0; i < omm_per_eye; ++i) {
                this->eye_plane.p.origin += this->position[i];
                this->eye_plane.p.normal += this->orientation[i];
            }
            this->eye_plane.p.origin /= omm_per_eye;
            this->eye_plane.p.normal.renormalize();
            this->eye_plane.ey = this->eye_plane.p.normal.cross (sm::vec<float>::uz());
            this->eye_plane.ey.renormalize();
            this->eye_plane.ex = this->eye_plane.ey.cross (this->eye_plane.p.normal); // guaranteed unit vector

            this->eye_plane.eye_to_oces.frombasis_inplace (this->eye_plane.ex, this->eye_plane.ey, this->eye_plane.p.normal);
            // That covers the rotation between the bases, but we also need a translation:
            this->eye_plane.eye_to_oces.pretranslate (this->eye_plane.p.origin);
        }

        // Find the maximum fields of view both horizontal and vertical
        void compute_fov_max()
        {
            auto horz_fov_r = sm::interval<float>::search_initialized();
            auto vert_fov_r = sm::interval<float>::search_initialized();

            const std::uint32_t omm_per_eye = this->get_omm_per_eye();

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

        void postprocess()
        {
            this->compute_fov_max();
            this->compute_neighbour_distance();
            this->compute_eye_plane();
            this->compute_projections_on_eye_plane();
            this->compute_central();
            this->ready = true;
        }

        // Make the mirror eye
        void construct_mirror_eye()
        {
            if (!this->mirrorplanes.empty()) {

                this->mirrors.clear();

                // Get size for one eye
                size_t sz = this->position.size();

                // Make space
                this->position.resize (2 * sz);
                this->orientation.resize (2 * sz);
                this->focal_offset.resize (2 * sz);
                this->diameter.resize (2 * sz);
                this->acceptance_angle.resize (2 * sz);

                sm::mat<float, 4> mirror = sm::mat<float, 4>::reflection (this->mirrorplanes[0].origin, this->mirrorplanes[0].normal);
                this->mirrors.push_back (mirror); // Saved for client code to use
                for (size_t i = 0; i < sz; ++i) {
                    // Mirror position and direction
                    sm::vec<float> mpos = (mirror * this->position[i]).less_one_dim();
                    sm::vec<float> mdir = (mirror * this->orientation[i]).less_one_dim();
                    this->position[sz + i] = mpos;
                    this->orientation[sz + i] = mdir;
                    // Focal offset, diameter and acceptance angle are simply copied
                    this->focal_offset[sz + i] = this->focal_offset[i];
                    this->diameter[sz + i] = this->diameter[i];
                    this->acceptance_angle[sz + i] = this->acceptance_angle[i];
                }
            }
        }

        void output_compound_ray_csv (const std::string& filename)
        {
            if (this->position.size() == this->orientation.size()
                && this->position.size() == this->focal_offset.size()
                && this->position.size() == this->acceptance_angle.size()) {

                std::ofstream fout (filename.c_str(), std::ios::out | std::ios::trunc);
                if (fout.is_open()) {
                    for (size_t i = 0; i < this->position.size(); ++i) {
                        fout << this->position[i].str_comma_separated (' ') << " "
                             << this->orientation[i].str_comma_separated (' ')
                             << " " << this->acceptance_angle[i]
                             << " " << std::abs(this->focal_offset[i])
                             << std::endl;
                    }
                    fout.close();
                }
            } else {
                std::cerr << "position, orientation, focal_offset and "
                          << "diameter/acceptance_angle should all have the same number of elements.\n";
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
}
