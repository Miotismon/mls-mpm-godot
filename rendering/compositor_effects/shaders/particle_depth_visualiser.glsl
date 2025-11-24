#[compute]
#version 450

layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;

layout(push_constant, std430) uniform Params {
    float particle_pos_tex_width;
    float particle_sphere_radius;
    vec2 render_size;

    float depth_divisor;
    // float padding;
    // float padding;
    // float padding;
};

layout(set = 0, binding = 0) uniform sampler2D fluid_depth_sampler;

layout(set = 1, binding = 0, rgba16f) restrict writeonly uniform image2D color_image;

layout(set = 2, binding = 0) uniform sampler2D fluid_linear_depth_sampler; 



layout(set = 3, binding = 0) uniform uniformBuffer {
    mat4 inv_proj;
};

void main()
{   

    ivec2 pixel_coord = ivec2(gl_GlobalInvocationID.xy);
    vec2 uv = (vec2(pixel_coord) + 0.5) / render_size; 
    vec2 texel_size = 1.0 / render_size;
	ivec2 size = ivec2(render_size);

	// Prevent reading/writing out of bounds.
	if (pixel_coord.x >= size.x || pixel_coord.y >= size.y) {
		return;
	}

	// // Read from our fluid color buffer.
    // vec2 depth_xy = imageLoad(fluid_color_image, pixel_coord).xy;

	// // Write back to our color buffer.
    // if (depth_xy != vec2(0.0))
    // {
    //     float depth = floor(depth_xy.x * 65536.0) / 65536.0 + floor(depth_xy.y * 65536.0) / (65536.0 * 65536.0);
    //     // imageStore(color_image, pixel_coord, vec4(vec3(depth), 1.0));
    //     float linear_depth = 1.0 / (depth * inv_proj[2].w + inv_proj[3].w);
    //     imageStore(color_image, pixel_coord, vec4(vec3(linear_depth / 100.0), 1.0));
        
    // }
	
    

    // Read from fluid linearized depth buffer
    float linear_depth = texture(fluid_linear_depth_sampler, uv, 0).r;

    if (linear_depth > 0.0 && linear_depth < 3990.0)
    {
        imageStore(color_image, pixel_coord, vec4(vec3(linear_depth / depth_divisor), 1.0));
    }

    // // Read from fluid depth buffer.
    // float depth = texelFetch(fluid_depth_sampler, pixel_coord, 0).r;

    // if (depth > 0.00001)
    // {
    //     // imageStore(color_image, pixel_coord, vec4(vec3(depth), 1.0));
    //     float linear_depth = 1.0 / (depth * inv_proj[2].w + inv_proj[3].w);
    //     imageStore(color_image, pixel_coord, vec4(vec3(linear_depth / 100.0), 1.0));
    // }

}

