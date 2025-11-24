shader_type spatial;
render_mode fog_disabled;


uniform float sphere_radius = 0.5;
uniform vec3 sphere_color : source_color = vec3(1.0);

uniform sampler2D particle_pos_tex;
uniform uint particle_pos_tex_width;

//varying vec3 particle_viewspace_pos;

void vertex() {
	
	// get particle position
	int idx = INSTANCE_ID;
	ivec2 coord = ivec2(idx % int(particle_pos_tex_width), idx / int(particle_pos_tex_width));
	vec3 particle_pos = texelFetch(particle_pos_tex, coord, 0).rgb;
	//VERTEX.xyz += particle_pos;
	VERTEX *= sphere_radius;
	//if (idx == 0)
	//{
		//particle_viewspace_pos = (VIEW_MATRIX * vec4(particle_pos, 1.0)).xyz;
	//}
	
	
	// example billboarding
	//MODELVIEW_MATRIX = VIEW_MATRIX * mat4(
			//MAIN_CAM_INV_VIEW_MATRIX[0],
			//MAIN_CAM_INV_VIEW_MATRIX[1],
			//MAIN_CAM_INV_VIEW_MATRIX[2],
			//MODEL_MATRIX[3]);
	//MODELVIEW_NORMAL_MATRIX = mat3(MODELVIEW_MATRIX);
	// custom billboarding
	MODELVIEW_MATRIX = VIEW_MATRIX * mat4(
		MAIN_CAM_INV_VIEW_MATRIX[0],
		MAIN_CAM_INV_VIEW_MATRIX[1],
		MAIN_CAM_INV_VIEW_MATRIX[2],
		MODEL_MATRIX[3] + vec4(particle_pos, 0.0) // add particle_pos
	);
	MODELVIEW_NORMAL_MATRIX = mat3(MODELVIEW_MATRIX);
	
	
}

void fragment() {
	
	// convert UV into [-1,1] range
	vec2 uv = UV * 2.0 - 1.0;
	float r2 = dot(uv, uv);

	// kill pixels outside the circle
	if (r2 > 1.0) {
	    discard;
	}
	
	// calculate normal in view space to fake spherical surface
	vec3 normal = normalize(vec3(uv.x, -uv.y, sqrt(1.0 - r2)));
	
	// calculate rotation matrix to rotate normals toward VIEW
	// Rodrigues' rotation formula
	vec3 v = vec3(-VIEW.y, VIEW.x, 0.0); // cross(vec3(0.0, 0.0, 1.0), VIEW)
	float s2 = v.x * v.x + v.y * v.y; // length(v) * length(v)
	float c = VIEW.z; // dot(vec3(0.0, 0.0, 1.0), VIEW);
	mat3 k_mat = mat3(vec3(0.0, v.z, -v.y), vec3(-v.z, 0.0, v.x), vec3(v.y, -v.x, 0.0));
	mat3 rot_mat = mat3(1.0) + k_mat + k_mat * k_mat * ((1.0 - c) / s2);
	
	NORMAL = rot_mat * normal; //NORMAL in fragment shader wants normal in view space
	
	
	// default depth calc
	//vec4 clip_pos = PROJECTION_MATRIX * vec4(VERTEX, 1.0);
	//clip_pos.xyz /= clip_pos.w;
	//DEPTH = clip_pos.z;
	// custom depth calc
	vec4 clip_pos = PROJECTION_MATRIX * vec4(VERTEX + normal.z * sphere_radius, 1.0); //VERTEX = interpolated view space pos of fragment
	//vec4 clip_pos = PROJECTION_MATRIX * vec4(particle_viewspace_pos + (normal) * sphere_radius, 1.0);
	DEPTH = clip_pos.z / clip_pos.w;
	
	// assign color
	ALBEDO = sphere_color;
	
}

void light() {
	vec3 new_diffuse = clamp(dot(NORMAL, LIGHT), 0.0, 1.0) * ATTENUATION * LIGHT_COLOR / PI;
	DIFFUSE_LIGHT += new_diffuse;
}
