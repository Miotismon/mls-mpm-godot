#[compute]
#version 450

layout(local_size_x = 8, local_size_y = 8, local_size_z = 1) in;

// Index of refraction for water
#define IOR 1.333

// Ratios of air and water IOR for refraction
// Air to water
#define ETA 1.0/IOR

// Fresnel at 0 degrees
#define F0 0.02

layout(set = 0, binding = 0) uniform sampler2D fluid_linear_depth_sampler; 

layout(set = 1, binding = 0) uniform sampler2D bg_depth_sampler;

layout(set = 2, binding = 0) uniform sampler2D src_color_sampler;

layout(set = 3, binding = 0, rgba16f) restrict writeonly uniform image2D dst_color_image;

layout(set = 4, binding = 0) uniform uniformBuffer {
    mat4 proj;

    mat4 inv_proj;
};

layout(set = 5, binding = 0) uniform samplerCube reflection_sampler_cube;

layout(push_constant, std430) uniform Params {
    vec2 render_size;
    float thickness;
    float optical_density;

    vec3 light_dir;
    float refraction_strength;

    vec3 diffuse_color;
    float specular_power;

    float fresnel_clamp;
    //float padding;
    //float padding;
    //float padding;
};


// float fresnelClamp;



vec3 get_view_pos_from_uv_and_depth(vec2 uv, float linear_depth)
{
    vec4 ndc;
    // ndc = vec4((uv * 2.0) - 1.0, depth, 1.0);
    ndc.xy = (uv * 2.0) - 1.0;
    ndc.z = proj[2].z + proj[3].z / linear_depth;
    ndc.w = 1.0;
    vec4 view_pos = inv_proj * ndc;
    return view_pos.xyz / view_pos.w;
}

vec3 get_view_pos_from_uv(vec2 uv) 
{
    float linear_depth = texture(fluid_linear_depth_sampler, uv, 0).r;
    return get_view_pos_from_uv_and_depth(uv, linear_depth);
}

void main()
{
    ivec2 pixel_coord = ivec2(gl_GlobalInvocationID.xy);
    vec2 uv = (vec2(pixel_coord) + vec2(0.5)) / render_size;
    vec2 texel_size = 1.0 / render_size;
    ivec2 size = ivec2(render_size);

    // Prevent reading/writing out of bounds.
    if(pixel_coord.x >= size.x || pixel_coord.y >= size.y) 
    {
        return;
    }

    float linear_depth = texture(fluid_linear_depth_sampler, uv, 0).r;

    if(linear_depth > 3990.0) 
    {
        return;
    }

    float bg_depth = texture(bg_depth_sampler, uv, 0).r;
    float linear_bg_depth = 1.0 / (bg_depth * inv_proj[2].w + inv_proj[3].w);

    if(linear_bg_depth < linear_depth) // if main render depth is infront of fluid depth
    {
        return;
    }
    

    vec4 bg_color = texture(src_color_sampler, uv, 0);

    // calculate view-space position of current pixel
    vec3 view_pos = get_view_pos_from_uv_and_depth(uv, linear_depth);

    // calculate differences 
    vec3 ddx = get_view_pos_from_uv(uv + vec2(texel_size.x, 0.)) - view_pos;
    vec3 ddx2 = view_pos - get_view_pos_from_uv(uv + vec2(-texel_size.x, 0.));
    if(abs(ddx.z) > abs(ddx2.z)) {
        ddx = ddx2;
    }

    vec3 ddy = get_view_pos_from_uv(uv + vec2(0., texel_size.y)) - view_pos;
    vec3 ddy2 = view_pos - get_view_pos_from_uv(uv + vec2(0., -texel_size.y));
    if(abs(ddy.z) > abs(ddy2.z)) {
        ddy = ddy2;
    }

    // calculate normal
    vec3 normal = normalize(cross(ddy, ddx));

    //imageStore(dst_color_image, pixel_coord, vec4(vec3(linear_depth / 100.0), 1.0));
    // imageStore(dst_color_image, pixel_coord, vec4(normal, 1.0));
    // return;
    
    // shading

    vec3 ray_dir = normalize(view_pos); // direction from camera position to view position

    // specular
    vec3  H = normalize(light_dir - ray_dir); // negative ray_dir because both vectors have to face away from the surface
    float specular = pow(max(0.0, dot(H, normal)), specular_power);


    // schlick's approximation
    float fresnel = clamp(F0 + (1.0 - F0) * pow(1.0 - dot(normal, -ray_dir), 5.0), 0., fresnel_clamp);

    // reflection of the cubemap
    vec3 reflection_dir = reflect(ray_dir, normal);
    vec3 reflection_color = (texture(reflection_sampler_cube, reflection_dir).rgb);

    // refraction of background render
    vec3 refraction_dir = refract(ray_dir, normal, ETA);

    vec4 transmitted = texture(src_color_sampler, vec2(uv + refraction_dir.xy * thickness * refraction_strength));

    vec3 transmittance = exp(-optical_density * (1.0 - diffuse_color) * thickness ); // Beer's law
   
    vec3 refraction_color = transmitted.rgb * transmittance;

    
    // combine refraction and reflection 
    vec3 final_color = mix(refraction_color, reflection_color, fresnel) + specular;


    //without reflections
    //vec3 final_color = refraction_color + specular;

    imageStore(dst_color_image, pixel_coord, vec4(final_color, transmitted.a));
}