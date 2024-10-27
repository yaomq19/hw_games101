#include <iostream>
#include <fstream>

#include <eigen3/Eigen/Eigen>
#include <lodepng.h>

#include "Triangle.hpp"
#include "rasterizer.hpp"
constexpr double MY_PI = 3.1415926;

Eigen::Matrix4f get_view_matrix(Eigen::Vector3f eye_pos)
{
    Eigen::Matrix4f view = Eigen::Matrix4f::Identity();

    Eigen::Matrix4f translate;
    translate << 1, 0, 0, -eye_pos[0], 0, 1, 0, -eye_pos[1], 0, 0, 1,
        -eye_pos[2], 0, 0, 0, 1;

    view = translate * view;

    return view;
}

Eigen::Matrix4f get_model_matrix(float rotation_angle)
{
     Eigen::Matrix4f model = Eigen::Matrix4f::Identity();
 
    // TODO: Implement this function
    // Create the model matrix for rotating the triangle around the Z axis.
    // Then return it.
    float angle = rotation_angle / 180 * MY_PI;
    model << std::cos(angle), -std::sin(angle), 0, 0,
        std::sin(angle), std::cos(angle), 0, 0,
        0, 0, 1, 0,
        0, 0, 0, 1;
 
    return model;
}

Eigen::Matrix4f get_projection_matrix(float eye_fov, float aspect_ratio,
                                      float zNear, float zFar)
{
    // Students will implement this function
    Eigen::Matrix4f projection = Eigen::Matrix4f::Identity();
    Eigen::Matrix4f Mperspective;
    Mperspective << zNear, 0, 0, 0,
        0, zNear, 0, 0,
        0, 0, zNear + zFar, -zNear * zFar,
        0, 0, 1, 0;
     
    /* 假定 eye_fov 是上下的角度 */
    /* zNear 需要取反，因为推导的矩阵是建立在 zNear ~ zFar 为负值的情况 */
    float half_height = std::tan(eye_fov / 2) * -zNear;
    float half_width = half_height * aspect_ratio;
 
    // 先平移后缩放，正交投影
    Eigen::Matrix4f Morth;
    Morth << 1 / half_width, 0, 0, 0,
        0, 1 / half_height, 0, 0,
        0, 0, 2 / (zNear - zFar), (zFar - zNear) / (zNear - zFar),
        0, 0, 0, 1;
 
    projection =   Morth * Mperspective;
    return projection;
}
std::vector<unsigned char> convert_vec_vec3f_to_vec_uchar(const std::vector<Eigen::Vector3f>& data) {
    std::vector<unsigned char> raw;
    raw.resize(data.size() * 3);
    for(int i = 0; i < data.size(); i++) {
        raw[i * 3] = data[i].x();
        raw[i * 3 + 1] = data[i].y();
        raw[i * 3 + 2] = data[i].z();
    }
    return raw;
}
int main(int argc, const char** argv)
{
    float angle = 190;

    rst::rasterizer r(700, 700);

    Eigen::Vector3f eye_pos = {0, 0, 5};

    std::vector<Eigen::Vector3f> pos{{2, 0, -2}, {0, 2, -2}, {-2, 0, -2}};

    std::vector<Eigen::Vector3i> ind{{0, 1, 2}};

    auto pos_id = r.load_positions(pos);
    auto ind_id = r.load_indices(ind);

    int key = 0;
    int frame_count = 0;

    r.clear(rst::Buffers::Color | rst::Buffers::Depth);

    r.set_model(get_model_matrix(angle));
    r.set_view(get_view_matrix(eye_pos));
    r.set_projection(get_projection_matrix(45, 1, 0.1, 50));

    r.draw(pos_id, ind_id, rst::Primitive::Triangle);

    std::vector<unsigned char> raw = convert_vec_vec3f_to_vec_uchar(r.frame_buffer());
    lodepng_encode_file("test.png", raw.data(), 700, 700, LCT_RGB, 8);                    
    return 0;
    
}
