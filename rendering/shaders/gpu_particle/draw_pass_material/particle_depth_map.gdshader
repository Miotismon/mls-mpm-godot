shader_type spatial;
render_mode fog_disabled, unshaded;

uniform float sphere_radius = 1.0;

void vertex() {

	// Don't render if not on the fluid layer
	if (CAMERA_VISIBLE_LAYERS != uint(2))
	{
		// This discards by setting the vertex to vec3(inf, inf, inf)
		VERTEX = vec3(1.0 / 0.0);
	}
	else
	{
		//VERTEX.xyz += particle_pos;
		VERTEX *= sphere_radius;
		//if (idx == 0)
		//{
			//particle_viewspace_pos = (VIEW_MATRIX * vec4(particle_pos, 1.0)).xyz;
		//}


		// example billboarding
		MODELVIEW_MATRIX = VIEW_MATRIX * mat4(
				MAIN_CAM_INV_VIEW_MATRIX[0],
				MAIN_CAM_INV_VIEW_MATRIX[1],
				MAIN_CAM_INV_VIEW_MATRIX[2],
				MODEL_MATRIX[3]);
		MODELVIEW_NORMAL_MATRIX = mat3(MODELVIEW_MATRIX);
		// custom billboarding
		//MODELVIEW_MATRIX = VIEW_MATRIX * mat4(
			//MAIN_CAM_INV_VIEW_MATRIX[0],
			//MAIN_CAM_INV_VIEW_MATRIX[1],
			//MAIN_CAM_INV_VIEW_MATRIX[2],
			//MODEL_MATRIX[3] + vec4(particle_pos, 0.0) // add particle_pos
		//);
		//MODELVIEW_NORMAL_MATRIX = mat3(MODELVIEW_MATRIX);
	}
}

void fragment() {

	// convert UV into [-1,1] range so 0,0 is center of circle
	vec2 uv = UV * 2.0 - 1.0;
	float r2 = dot(uv, uv);

	// kill pixels outside the circle
	if (r2 > 1.0) {
		DEPTH = 0.0;
	    discard;
	}

	// calculate normal in view space to fake spherical surface
	vec3 normal = normalize(vec3(uv.x, -uv.y, sqrt(1.0 - r2)));

	// default depth calc
	//vec4 clip_pos = PROJECTION_MATRIX * vec4(VERTEX, 1.0);
	//clip_pos.xyz /= clip_pos.w;
	//DEPTH = clip_pos.z;

	// custom depth calc
	//vec4 clip_pos = PROJECTION_MATRIX * vec4(particle_viewspace_pos + (normal) * sphere_radius, 1.0);
	vec4 clip_pos = PROJECTION_MATRIX * vec4(VERTEX + normal.z * sphere_radius, 1.0); //VERTEX = interpolated view space pos of fragment
	float depth = clip_pos.z / clip_pos.w;
	DEPTH = depth;

	//draw depth into visual layer 2
	ALBEDO.x = depth;
	ALBEDO.y = fract(depth * 65536.0);
	ALBEDO.z = 0.0; // Could add other info in this channel if desired...

	// debug see depth
	//vec3 ndc = vec3(SCREEN_UV * 2.0 - 1.0, DEPTH);
	//vec4 view = INV_PROJECTION_MATRIX * vec4(ndc, 1.0); // this is for calculating view space pos from depth, here only view.z is needed so this can be simplified (https://www.youtube.com/watch?v=NCXr8zrT5zs)
	//view.xyz /= view.w;
	//float linear_depth = -view.z; // depth in view space in meters, range[0, max view dist], 0 is touching camera
	//ALBEDO = vec3(linear_depth) / 100.0; // visualising depth between 0 and 100 meters as albedo 0 to 1
}
