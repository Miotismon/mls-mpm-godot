shader_type spatial;
render_mode unshaded, fog_disabled;

uniform sampler2D screen_texture : source_color, hint_screen_texture;
uniform sampler2D depth_texture : hint_depth_texture;
//uniform sampler2D normal_rough_texture : hint_normal_roughness_texture;

uniform int mode : hint_enum("Depth", "Color") = 0;

// x or y
uniform vec2 blur_dir = vec2(1.0, 0.0);

//uniform float particle_radius = 0.5;
//uniform int filter_size = 10;
//uniform int max_filter_size = 100;
//const float thresholdRatio = 10.0;



void vertex() {
	POSITION = vec4(VERTEX.xy, 1.0, 1.0);

	//projectedParticleConstant = (blurFilterSize * sphere_radius * 2.0 * 0.05 * (VIEWPORT_SIZE.y / 2.0)) / tan((5.0*PI / 12.0) / 2.0);
}

//float compute_weight2D(vec2 r, float two_sigma2)
//{
    //return exp(-dot(r, r) / two_sigma2);
//}
//
//float compute_weight1D(float r, float two_sigma2)
//{
    //return exp(-r * r / two_sigma2);
//}

void fragment() {
	vec2 pixel_size = 1.0 / VIEWPORT_SIZE;
	
	float depth = texture(depth_texture, SCREEN_UV).r;
	float linear_depth = 1.0 / (depth * INV_PROJECTION_MATRIX[2].w + INV_PROJECTION_MATRIX[3].w);
	
	if (depth <= 0.0001) // 0.0 == far plane, also direct floating-point comparison (depth == 0.0) is bad
	{
		//ALBEDO = vec3(depth);
		discard;
	}
	else if (mode == 0) // Depth mode
	{
		
		int filter_size = 50;
		float sigma_space = float(filter_size) / 3.0;
		float two_sigma_space2 = 2.0 * sigma_space * sigma_space;

		float sigma_range = 2.0;
		float two_sigma_range2 = 2.0 * sigma_range * sigma_range;

		float sum = 0.0;
		float wsum = 0.0;

		for (int x = -filter_size; x <= filter_size; ++x) 
		{
			float sample = textureLod(depth_texture, SCREEN_UV + (float(x) * pixel_size * blur_dir), 0.).r;
			float linear_sample = 1.0 / (sample * INV_PROJECTION_MATRIX[2].w + INV_PROJECTION_MATRIX[3].w);
			
			float r = float(x * x);
			float w = exp(-r / two_sigma_space2);

			float rDepth = linear_sample - linear_depth;
			float wd = exp(-rDepth * rDepth / two_sigma_range2);

			sum += linear_sample * w * wd;
			wsum += w * wd;
		}
		sum /= wsum;
		
		
		
		float test_depth = texture(depth_texture, SCREEN_UV).r;
		float test_linear_depth = 1.0 / (test_depth * INV_PROJECTION_MATRIX[2].w + INV_PROJECTION_MATRIX[3].w);
		//float test_linear_depth_div_100 = test_linear_depth / 100.0;
		
		//DEPTH = clamp(sum, 0.0, 1.0);
		//float linear_sum = 1.0 / (sum * INV_PROJECTION_MATRIX[2].w + INV_PROJECTION_MATRIX[3].w);
		//float clip_sum = (PROJECTION_MATRIX[2].z * sum + PROJECTION_MATRIX[3].z) / -sum; 
		//DEPTH = clip_sum;
		ALBEDO = vec3(sum / 100.0);
		//ALBEDO = vec3(test_linear_depth / 100.0);
		//ALBEDO = texture(screen_texture, SCREEN_UV).rgb;
	}
	else if (mode == 1) // color mode
	{
		vec3 source = textureLod(screen_texture, SCREEN_UV, 0.).rgb;

		int filter_size = 30;
		float sigma = float(filter_size) / 3.0;
		float two_sigma2 = 2.0 * sigma * sigma;

		float sigma_range = 1.0;
		float two_sigma_range2 = 2.0 * sigma_range * sigma_range;

		vec3 sum = vec3(0.0);
		vec3 wsum = vec3(0.0);

		for (int x = -filter_size; x <= filter_size; ++x) 
		{
			vec3 sample = textureLod(screen_texture, SCREEN_UV + (float(x) * pixel_size * blur_dir), 0.).rgb;

			float r = float(x * x);
			float w = exp(-r / two_sigma2);

			vec3 rDepth = sample.r - source;
			vec3 wd = exp(-rDepth * rDepth / two_sigma_range2);

			sum += sample * w * wd;
			wsum += w * wd;
		}
		sum /= wsum;
		
		
		ALBEDO = sum;
		//ALBEDO = texture(screen_texture, SCREEN_UV).rgb;
		
	}

}

//void light() {
//	// Called for every pixel for every light affecting the material.
//	// Uncomment to replace the default light processing function with this one.
//}
