
extends CompositorEffect
class_name ScreenSpaceFluidRendering

enum RenderType {DEFAULT, VELOCITY_SPHERES}

@export var render_type: RenderType = RenderType.DEFAULT

@export var particle_sphere_radius: float = 1.0

@export_group("Depth Blur")
@export var depth_blur_enabled: bool = true
@export var blur_depth_scale: float = 10.0
@export var max_filter_size: int = 100
@export var blur_filter_size: float = 7.0

@export_group("Fluid Render")
@export var diffuse_color: Color = Color(0.085, 0.6375, 0.765)
@export_range(0.0, 3.0) var minimum_thickness: float = 1.0
@export var optical_density: float = 2.0
##Generally between 0 and 0.3
@export_range(0.0, 0.5) var refraction_strength: float = 0.1
##Increase the value to make the specular effect more concentrated
@export var specular_power: float = 250.0
@export var fresnel_clamp: float = 1.0

@export_group("Debug Render")
@export var debug_draw_depth: bool = false
@export var depth_divisor: float = 100.0

var light_world_dir: Vector3 = Vector3(0.0, -1.0, 0.0)


var fluid_color_image: RID
var fluid_depth_image: RID

var rd: RenderingDevice

var nearest_sampler: RID
var linear_sampler: RID

var linearize_depth_shader: RID
var linearize_depth_pipeline: RID 

var bilateral_blur_shader: RID
var bilateral_blur_pipeline: RID

var fluid_render_fixed_depth_shader: RID
var fluid_render_fixed_depth_pipeline: RID

var copy_texture_shader: RID
var copy_texture_pipeline: RID

var particle_depth_visualiser_shader: RID
var particle_depth_visualiser_pipeline: RID

var fluid_render_velocity_spheres_shader: RID
var fluid_render_velocity_spheres_pipeline: RID

var context: StringName = "ScreenSpaceFluidRendering"
var fluid_linear_depth_texture: StringName = "fluid_linear_depth_texture"
var temp_fluid_linear_depth_texture: StringName = "temp_fluid_linear_depth_texture"
var temp_color_texture: StringName = "temp_color_texture"
#var reflections_cubemap_texture: StringName = "reflections_cubemap_texture"
#var reflections_cubemap_1: CompressedTexture2D = preload("res://assets/cubemap/cubemap_template1.png")
#var reflections_cubemap_2: CompressedTexture2D = preload("res://assets/cubemap/cubemap_template2.png")
#var reflections_cubemap_3: CompressedTexture2D = preload("res://assets/cubemap/cubemap_template3.png")
#var reflections_cubemap_4: CompressedTexture2D = preload("res://assets/cubemap/cubemap_template4.png")
#var reflections_cubemap_5: CompressedTexture2D = preload("res://assets/cubemap/cubemap_template5.png")
#var reflections_cubemap_6: CompressedTexture2D = preload("res://assets/cubemap/cubemap_template6.png")
var reflections_cubemap_1: CompressedTexture2D = preload("res://assets/cubemap/Cubemap_Sky_04-1.png")
var reflections_cubemap_2: CompressedTexture2D = preload("res://assets/cubemap/Cubemap_Sky_04-2.png")
var reflections_cubemap_3: CompressedTexture2D = preload("res://assets/cubemap/Cubemap_Sky_04-3.png")
var reflections_cubemap_4: CompressedTexture2D = preload("res://assets/cubemap/Cubemap_Sky_04-4.png")
var reflections_cubemap_5: CompressedTexture2D = preload("res://assets/cubemap/Cubemap_Sky_04-5.png")
var reflections_cubemap_6: CompressedTexture2D = preload("res://assets/cubemap/Cubemap_Sky_04-6.png")
var reflections_cubemap_1_image: Image = reflections_cubemap_1.get_image()
var reflections_cubemap_2_image: Image = reflections_cubemap_2.get_image()
var reflections_cubemap_3_image: Image = reflections_cubemap_3.get_image()
var reflections_cubemap_4_image: Image = reflections_cubemap_4.get_image()
var reflections_cubemap_5_image: Image = reflections_cubemap_5.get_image()
var reflections_cubemap_6_image: Image = reflections_cubemap_6.get_image()
var cubemap_texture: RID

var inv_projection_matrix_uniform_buffer : RID
var inv_projection_projection_matrix_uniform_buffer : RID

func _init() -> void:
	effect_callback_type = EFFECT_CALLBACK_TYPE_POST_TRANSPARENT
	rd = RenderingServer.get_rendering_device()
	RenderingServer.call_on_render_thread(initialize_compute)

# System notifications, we want to react on the notification that alerts us we are about to be destroyed.
func _notification(what : int) -> void:
	if what == NOTIFICATION_PREDELETE:
		
		if nearest_sampler.is_valid():
			rd.free_rid(nearest_sampler)
		if linear_sampler.is_valid():
			rd.free_rid(linear_sampler)
		
		# Freeing our shader will also free any dependents such as the pipeline!
		if linearize_depth_shader.is_valid():
			rd.free_rid(linearize_depth_shader)
		if bilateral_blur_shader.is_valid():
			rd.free_rid(bilateral_blur_shader)
		if fluid_render_fixed_depth_shader.is_valid():
			rd.free_rid(fluid_render_fixed_depth_shader)
		if copy_texture_shader.is_valid():
			rd.free_rid(copy_texture_shader)
		if particle_depth_visualiser_shader.is_valid():
			rd.free_rid(particle_depth_visualiser_shader)
		if fluid_render_velocity_spheres_shader.is_valid():
			rd.free_rid(fluid_render_velocity_spheres_shader)


#region Code in this region runs on the rendering thread.
# Compile shaders at initialization.
func initialize_compute() -> void:
	rd = RenderingServer.get_rendering_device()
	if not rd:
		return
	
	var sampler_state : RDSamplerState = RDSamplerState.new()
	sampler_state.min_filter = RenderingDevice.SAMPLER_FILTER_NEAREST
	sampler_state.mag_filter = RenderingDevice.SAMPLER_FILTER_NEAREST
	nearest_sampler = rd.sampler_create(sampler_state)
	
	sampler_state = RDSamplerState.new()
	sampler_state.min_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	sampler_state.mag_filter = RenderingDevice.SAMPLER_FILTER_LINEAR
	linear_sampler = rd.sampler_create(sampler_state)
	
	# Compile shaders.
	var shader_file : RDShaderFile = load("res://rendering/compositor_effects/shaders/bilateral_blur_directional.glsl")
	var shader_spirv : RDShaderSPIRV = shader_file.get_spirv()
	bilateral_blur_shader = rd.shader_create_from_spirv(shader_spirv)
	if bilateral_blur_shader.is_valid():
		bilateral_blur_pipeline = rd.compute_pipeline_create(bilateral_blur_shader)
	
	shader_file = load("res://rendering/compositor_effects/shaders/linearize_fluid_depth.glsl")
	shader_spirv = shader_file.get_spirv()
	linearize_depth_shader = rd.shader_create_from_spirv(shader_spirv)
	if linearize_depth_shader.is_valid():
		linearize_depth_pipeline = rd.compute_pipeline_create(linearize_depth_shader)
	
	shader_file = load("res://rendering/compositor_effects/shaders/fluid_render_fixed_depth.glsl")
	shader_spirv = shader_file.get_spirv()
	fluid_render_fixed_depth_shader = rd.shader_create_from_spirv(shader_spirv)
	if fluid_render_fixed_depth_shader.is_valid():
		fluid_render_fixed_depth_pipeline = rd.compute_pipeline_create(fluid_render_fixed_depth_shader)
	
	shader_file = load("res://rendering/compositor_effects/shaders/copy_texture.glsl")
	shader_spirv = shader_file.get_spirv()
	copy_texture_shader = rd.shader_create_from_spirv(shader_spirv)
	if copy_texture_shader.is_valid():
		copy_texture_pipeline = rd.compute_pipeline_create(copy_texture_shader)
	
	shader_file = load("res://rendering/compositor_effects/shaders/particle_depth_visualiser.glsl")
	shader_spirv = shader_file.get_spirv()
	particle_depth_visualiser_shader = rd.shader_create_from_spirv(shader_spirv)
	if particle_depth_visualiser_shader.is_valid():
		particle_depth_visualiser_pipeline = rd.compute_pipeline_create(particle_depth_visualiser_shader)
	
	shader_file = load("res://rendering/compositor_effects/shaders/fluid_render_velocity_spheres.glsl")
	shader_spirv = shader_file.get_spirv()
	fluid_render_velocity_spheres_shader = rd.shader_create_from_spirv(shader_spirv)
	if fluid_render_velocity_spheres_shader.is_valid():
		fluid_render_velocity_spheres_pipeline = rd.compute_pipeline_create(fluid_render_velocity_spheres_shader)
	
	var usage_bits : int = RenderingDevice.TEXTURE_USAGE_SAMPLING_BIT | RenderingDevice.TEXTURE_USAGE_STORAGE_BIT | RenderingDevice.TEXTURE_USAGE_CAN_COPY_TO_BIT | RenderingDevice.TEXTURE_USAGE_CAN_COPY_FROM_BIT
	var cubemap_format: RDTextureFormat = RDTextureFormat.new()
	cubemap_format.array_layers = 6
	cubemap_format.format = RenderingDevice.DATA_FORMAT_R8G8B8A8_UNORM
	cubemap_format.texture_type = RenderingDevice.TEXTURE_TYPE_CUBE
	cubemap_format.usage_bits = usage_bits
	cubemap_format.height = 512
	cubemap_format.width = 512
	var cubemap_view: RDTextureView = RDTextureView.new()
	var cubemap_data: Array[PackedByteArray] = [reflections_cubemap_1_image.get_data(), reflections_cubemap_2_image.get_data(), reflections_cubemap_3_image.get_data(), reflections_cubemap_4_image.get_data(), reflections_cubemap_5_image.get_data(), reflections_cubemap_6_image.get_data()]
	#render_scene_buffers.create_texture_from_format(context, reflections_cubemap_texture, cube_map_format, cube_map_view, false)
	cubemap_texture = rd.texture_create(cubemap_format, cubemap_view, cubemap_data)
	print("reflections_cubemap_1_image format: ", reflections_cubemap_1_image.get_format())
	
	inv_projection_matrix_uniform_buffer = rd.uniform_buffer_create(64, [])
	inv_projection_projection_matrix_uniform_buffer = rd.uniform_buffer_create(128, [])

func validate_pipelines():
	return particle_depth_visualiser_pipeline.is_valid()

# Called by the rendering thread every frame.
func _render_callback(p_effect_callback_type: EffectCallbackType, p_render_data: RenderData) -> void:
	if rd and p_effect_callback_type == EFFECT_CALLBACK_TYPE_POST_TRANSPARENT and validate_pipelines():
		# Get our render scene buffers object, this gives us access to our render buffers.
		# Note that implementation differs per renderer hence the need for the cast.
		var render_scene_buffers : RenderSceneBuffersRD = p_render_data.get_render_scene_buffers()
		var render_scene_data : RenderSceneDataRD = p_render_data.get_render_scene_data()
		if render_scene_buffers and render_scene_data:
			# Get our render size, this is the 3D render resolution!
			var render_size : Vector2i = render_scene_buffers.get_internal_size()
			if render_size.x == 0 and render_size.y == 0:
				return

			# If we have buffers for this viewport, check if they are the right size
			if render_scene_buffers.has_texture(context, fluid_linear_depth_texture):
				var tf : RDTextureFormat = render_scene_buffers.get_texture_format(context, fluid_linear_depth_texture)
				if tf.width != render_size.x or tf.height != render_size.y:
					# This will clear all textures for this viewport under this context
					render_scene_buffers.clear_context(context)

			if !render_scene_buffers.has_texture(context, fluid_linear_depth_texture):
				var usage_bits : int = RenderingDevice.TEXTURE_USAGE_SAMPLING_BIT | RenderingDevice.TEXTURE_USAGE_STORAGE_BIT | RenderingDevice.TEXTURE_USAGE_CAN_COPY_TO_BIT | RenderingDevice.TEXTURE_USAGE_CAN_COPY_FROM_BIT
				render_scene_buffers.create_texture(context, fluid_linear_depth_texture, RenderingDevice.DATA_FORMAT_R32_SFLOAT, usage_bits, RenderingDevice.TEXTURE_SAMPLES_1, render_size, 1, 1, true, false)
				render_scene_buffers.create_texture(context, temp_fluid_linear_depth_texture, RenderingDevice.DATA_FORMAT_R32_SFLOAT, usage_bits, RenderingDevice.TEXTURE_SAMPLES_1, render_size, 1, 1, true, false)
				render_scene_buffers.create_texture(context, temp_color_texture, RenderingDevice.DATA_FORMAT_R16G16B16A16_SFLOAT, usage_bits, RenderingDevice.TEXTURE_SAMPLES_1, render_size, 1, 1, true, false)
				#var cubemap_format: RDTextureFormat = RDTextureFormat.new()
				#cubemap_format.array_layers = 6
				#cubemap_format.format = RenderingDevice.DATA_FORMAT_R16G16B16A16_SFLOAT
				#cubemap_format.texture_type = RenderingDevice.TEXTURE_TYPE_CUBE
				#cubemap_format.usage_bits = usage_bits
				#cubemap_format.height = 512
				#cubemap_format.width = 512
				#var cubemap_view: RDTextureView = RDTextureView.new()
				#var cubemap_data: Array[PackedByteArray] = [reflections_cubemap_1_image.get_data(), PackedByteArray(), PackedByteArray(), PackedByteArray(), PackedByteArray(), PackedByteArray()]
				##render_scene_buffers.create_texture_from_format(context, reflections_cubemap_texture, cube_map_format, cube_map_view, false)
				#var cubemap_texture = rd.texture_create(cubemap_format, cubemap_view, cubemap_data)

			# Loop through views just in case we're doing stereo rendering. No extra cost if this is mono.
			var view_count : int = render_scene_buffers.get_view_count()
			for view in view_count:
				# Get necessary images
				var color_image : RID = render_scene_buffers.get_color_layer(view)
				var depth_image : RID = render_scene_buffers.get_depth_layer(view)
				
				if !rd.texture_is_valid(fluid_color_image):
					print("getting fluid_color_image ref")
					fluid_color_image = Global.fluid_color_image
					if fluid_color_image.is_valid():
						print("Fluid color format:", rd.texture_get_format(fluid_color_image).format)
				
				if !rd.texture_is_valid(fluid_depth_image):
					print("getting fluid_depth_image ref")
					fluid_depth_image = Global.fluid_depth_image
					if fluid_depth_image.is_valid():
						print("Fluid depth format:", rd.texture_get_format(fluid_depth_image).format)
				
				var fluid_linear_depth_image: RID = render_scene_buffers.get_texture_slice(context, fluid_linear_depth_texture, view, 0, 1, 1)
				var temp_fluid_linear_depth_image: RID = render_scene_buffers.get_texture_slice(context, temp_fluid_linear_depth_texture, view, 0, 1, 1)
				var temp_color_image: RID = render_scene_buffers.get_texture_slice(context, temp_color_texture, view, 0, 1, 1)
				#var reflections_cubemap_image: RID = render_scene_buffers.get_texture(context, reflections_cubemap_texture)
				
				# clear custom textures
				rd.texture_clear(fluid_linear_depth_image, Color(), 0, 1, view, 1)
				rd.texture_clear(temp_fluid_linear_depth_image, Color(), 0, 1, view, 1)
				rd.texture_clear(temp_color_image, Color(0.0, 0.0, 0.0, 0.0), 0, 1, view, 1)
				
				#rd.texture_clear(reflections_cubemap_image, Color(1.0, 0.0, 0.0, 1.0), 0, 1, 0, 6)
				
				# Get some rendering info
				var cam_tr : Transform3D = render_scene_data.get_cam_transform()
				var view_tr : Transform3D = cam_tr.affine_inverse()
				var proj : Projection = render_scene_data.get_view_projection(view)
				var inv_proj : Projection = proj.inverse()
				
				var view_mat = [
					view_tr.basis.x.x, view_tr.basis.x.y, view_tr.basis.x.z, 0.0,
					view_tr.basis.y.x, view_tr.basis.y.y, view_tr.basis.y.z, 0.0,
					view_tr.basis.z.x, view_tr.basis.z.y, view_tr.basis.z.z, 0.0,
					view_tr.origin.x, view_tr.origin.y, view_tr.origin.z, 1.0,
				]
				
				var projection_mat = [
					proj.x.x, proj.x.y, proj.x.z, proj.x.w,
					proj.y.x, proj.y.y, proj.y.z, proj.y.w, 
					proj.z.x, proj.z.y, proj.z.z, proj.z.w,
					proj.w.x, proj.w.y, proj.w.z, proj.w.w,
				]
				
				var inv_projection_mat = [
					inv_proj.x.x, inv_proj.x.y, inv_proj.x.z, inv_proj.x.w,
					inv_proj.y.x, inv_proj.y.y, inv_proj.y.z, inv_proj.y.w, 
					inv_proj.z.x, inv_proj.z.y, inv_proj.z.z, inv_proj.z.w,
					inv_proj.w.x, inv_proj.w.y, inv_proj.w.z, inv_proj.w.w,
				]
				
				var _view_mat_array : PackedByteArray = PackedFloat32Array(view_mat).to_byte_array()
				var projection_mat_array : PackedByteArray = PackedFloat32Array(projection_mat).to_byte_array()
				var inv_projection_mat_array : PackedByteArray = PackedFloat32Array(inv_projection_mat).to_byte_array()
				
				var pb = PackedByteArray()
				pb.append_array(inv_projection_mat_array)
				rd.buffer_update(inv_projection_matrix_uniform_buffer, 0, 64, inv_projection_mat_array)
				
				pb = PackedByteArray()
				pb.append_array(projection_mat_array)
				pb.append_array(inv_projection_mat_array)
				rd.buffer_update(inv_projection_projection_matrix_uniform_buffer, 0, 128, pb)
				
				
				# start fetching/creating uniforms and push constants for each shader
				
				if render_type == RenderType.DEFAULT:
					
					# Step 0: linearize depth
					var uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(linear_sampler)
					uniform.add_id(fluid_depth_image)
					var fluid_depth_uniform_set := UniformSetCacheRD.get_cache(linearize_depth_shader, 0, [uniform])
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
					uniform.binding = 0
					uniform.add_id(fluid_linear_depth_image)
					var fluid_linear_depth_uniform_set := UniformSetCacheRD.get_cache(linearize_depth_shader, 1, [uniform])
					
					#uniform = RDUniform.new()
					#uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
					#uniform.binding = 0
					#uniform.add_id(temp_fluid_linear_depth_image)
					#var temp_fluid_linear_depth_uniform_set := UniformSetCacheRD.get_cache(linearize_depth_shader, 2, [uniform])
					
					
					var matrices_uniform : RDUniform = RDUniform.new()
					matrices_uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_UNIFORM_BUFFER
					matrices_uniform.binding = 0
					matrices_uniform.add_id(inv_projection_matrix_uniform_buffer)
					var matrices_uniform_set: RID = UniformSetCacheRD.get_cache(linearize_depth_shader, 3, [matrices_uniform])
					
					var push_constant : PackedFloat32Array = PackedFloat32Array()
					push_constant.push_back(render_size.x)
					push_constant.push_back(render_size.y)
					push_constant.push_back(0.0)
					push_constant.push_back(0.0)
					
					@warning_ignore("integer_division")
					var x_groups := (render_size.x - 1) / 8 + 1
					@warning_ignore("integer_division")
					var y_groups := (render_size.y - 1) / 8 + 1
					var z_groups := 1
					
					var compute_list := rd.compute_list_begin()
					rd.compute_list_bind_compute_pipeline(compute_list, linearize_depth_pipeline)
					rd.compute_list_bind_uniform_set(compute_list, fluid_depth_uniform_set, 0)
					rd.compute_list_bind_uniform_set(compute_list, fluid_linear_depth_uniform_set, 1)
					#rd.compute_list_bind_uniform_set(compute_list, temp_fluid_linear_depth_uniform_set, 2)
					rd.compute_list_bind_uniform_set(compute_list, matrices_uniform_set, 3)
					rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
					rd.compute_list_dispatch(compute_list, x_groups, y_groups, z_groups)
					rd.compute_list_end()
					
					rd.texture_copy(fluid_linear_depth_image, temp_fluid_linear_depth_image, Vector3.ZERO, Vector3.ZERO, Vector3(render_size.x, render_size.y, 0), 0, 0, 0, 0)
					
					# Step 1: bilateral blur
					if depth_blur_enabled:
						
						# Step 1.1: blur in x dir
						uniform = RDUniform.new()
						uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
						uniform.binding = 0
						uniform.add_id(linear_sampler)
						uniform.add_id(fluid_linear_depth_image)
						var src_fluid_linear_depth_uniform_set = UniformSetCacheRD.get_cache(bilateral_blur_shader, 0, [uniform])
						
						uniform = RDUniform.new()
						uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
						uniform.binding = 0
						uniform.add_id(temp_fluid_linear_depth_image)
						var dst_fluid_linear_depth_uniform_set = UniformSetCacheRD.get_cache(bilateral_blur_shader, 1, [uniform])
						
						# depth_threshold = (paritcle_diameter / 2) * blur_depth_scale;
						var depth_threshold : float = particle_sphere_radius * blur_depth_scale
						# max_filter_size from export
						# projected_particle_constant = (blur_filter_size * paritcle_diameter * 0.05 * (render_size.y / 2)) / tan(fov / 2)
						var projected_particle_constant = (blur_filter_size * particle_sphere_radius * 0.1 * (render_size.y / 2.0)) / tan((75 * PI) / 180 / 2.0)
						
						push_constant = PackedFloat32Array()
						push_constant.push_back(1.0)
						push_constant.push_back(0.0)
						push_constant.push_back(render_size.x)
						push_constant.push_back(render_size.y)
						push_constant.push_back(depth_threshold)
						push_constant.push_back(max_filter_size)
						push_constant.push_back(projected_particle_constant)
						push_constant.push_back(0.0)
						
						compute_list = rd.compute_list_begin()
						rd.compute_list_bind_compute_pipeline(compute_list, bilateral_blur_pipeline)
						rd.compute_list_bind_uniform_set(compute_list, src_fluid_linear_depth_uniform_set, 0)
						rd.compute_list_bind_uniform_set(compute_list, dst_fluid_linear_depth_uniform_set, 1)
						rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
						rd.compute_list_dispatch(compute_list, x_groups, y_groups, z_groups)
						rd.compute_list_end()
						
						
						# Step 1.2: blur in y dir
						uniform = RDUniform.new()
						uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
						uniform.binding = 0
						uniform.add_id(linear_sampler)
						uniform.add_id(temp_fluid_linear_depth_image)
						src_fluid_linear_depth_uniform_set = UniformSetCacheRD.get_cache(bilateral_blur_shader, 0, [uniform])
						
						uniform = RDUniform.new()
						uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
						uniform.binding = 0
						uniform.add_id(fluid_linear_depth_image)
						dst_fluid_linear_depth_uniform_set = UniformSetCacheRD.get_cache(bilateral_blur_shader, 1, [uniform])
						
						push_constant = PackedFloat32Array()
						push_constant.push_back(0.0)
						push_constant.push_back(1.0)
						push_constant.push_back(render_size.x)
						push_constant.push_back(render_size.y)
						push_constant.push_back(depth_threshold)
						push_constant.push_back(max_filter_size)
						push_constant.push_back(projected_particle_constant)
						push_constant.push_back(0.0)
						
						compute_list = rd.compute_list_begin()
						rd.compute_list_bind_compute_pipeline(compute_list, bilateral_blur_pipeline)
						rd.compute_list_bind_uniform_set(compute_list, src_fluid_linear_depth_uniform_set, 0)
						rd.compute_list_bind_uniform_set(compute_list, dst_fluid_linear_depth_uniform_set, 1)
						rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
						rd.compute_list_dispatch(compute_list, x_groups, y_groups, z_groups)
						rd.compute_list_end()
					
					# Step 2: get fluid depth (didn't have time for that)
					
					# Step 3: render out fluid into temp color texture
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(linear_sampler)
					uniform.add_id(fluid_linear_depth_image)
					fluid_linear_depth_uniform_set = UniformSetCacheRD.get_cache(fluid_render_fixed_depth_shader, 0, [uniform])
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(linear_sampler)
					uniform.add_id(depth_image)
					var bg_depth_uniform_set := UniformSetCacheRD.get_cache(fluid_render_fixed_depth_shader, 1, [uniform])
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(linear_sampler)
					uniform.add_id(color_image)
					var color_uniform_set := UniformSetCacheRD.get_cache(fluid_render_fixed_depth_shader, 2, [uniform])
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
					uniform.binding = 0
					uniform.add_id(temp_color_image)
					var temp_color_uniform_set := UniformSetCacheRD.get_cache(fluid_render_fixed_depth_shader, 3, [uniform])
					
					
					matrices_uniform = RDUniform.new()
					matrices_uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_UNIFORM_BUFFER
					matrices_uniform.binding = 0
					matrices_uniform.add_id(inv_projection_projection_matrix_uniform_buffer)
					matrices_uniform_set = UniformSetCacheRD.get_cache(fluid_render_fixed_depth_shader, 4, [matrices_uniform])
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(linear_sampler)
					uniform.add_id(cubemap_texture)
					var reflections_cubemap_uniform_set := UniformSetCacheRD.get_cache(fluid_render_fixed_depth_shader, 5, [uniform])
					
					
					var light_view_dir = view_tr.basis * light_world_dir
					light_view_dir = -light_view_dir.normalized() # has to point toward light source
					
					
					
					push_constant = PackedFloat32Array()
					push_constant.append(render_size.x)
					push_constant.append(render_size.y)
					push_constant.append(minimum_thickness) # thickness
					push_constant.append(optical_density) # optical_density
					push_constant.append(light_view_dir.x) # light_view_dir 
					push_constant.append(light_view_dir.y)
					push_constant.append(light_view_dir.z)
					push_constant.append(refraction_strength) # refraction_strength 
					push_constant.append(diffuse_color.r) # diffuse_color
					push_constant.append(diffuse_color.g)
					push_constant.append(diffuse_color.b)
					push_constant.append(specular_power) # specular_power
					push_constant.append(fresnel_clamp)
					push_constant.append(0.0)
					push_constant.append(0.0)
					push_constant.append(0.0)
					
					compute_list = rd.compute_list_begin()
					rd.compute_list_bind_compute_pipeline(compute_list, fluid_render_fixed_depth_pipeline)
					rd.compute_list_bind_uniform_set(compute_list, fluid_linear_depth_uniform_set, 0)
					rd.compute_list_bind_uniform_set(compute_list, bg_depth_uniform_set, 1)
					rd.compute_list_bind_uniform_set(compute_list, color_uniform_set, 2)
					rd.compute_list_bind_uniform_set(compute_list, temp_color_uniform_set, 3)
					rd.compute_list_bind_uniform_set(compute_list, matrices_uniform_set, 4)
					rd.compute_list_bind_uniform_set(compute_list, reflections_cubemap_uniform_set, 5)
					rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
					rd.compute_list_dispatch(compute_list, x_groups, y_groups, z_groups)
					rd.compute_list_end()
					
					# Step 4: copy temp color to real color
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(nearest_sampler)
					uniform.add_id(temp_color_image)
					temp_color_uniform_set = UniformSetCacheRD.get_cache(copy_texture_shader, 0, [uniform])
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
					uniform.binding = 0
					uniform.add_id(color_image)
					color_uniform_set = UniformSetCacheRD.get_cache(copy_texture_shader, 1, [uniform])
					
					push_constant = PackedFloat32Array()
					push_constant.append(render_size.x)
					push_constant.append(render_size.y)
					push_constant.append(0.000001) # color_threshold
					push_constant.append(0.0)
					
					compute_list = rd.compute_list_begin()
					rd.compute_list_bind_compute_pipeline(compute_list, copy_texture_pipeline)
					rd.compute_list_bind_uniform_set(compute_list, temp_color_uniform_set, 0)
					rd.compute_list_bind_uniform_set(compute_list, color_uniform_set, 1)
					rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
					rd.compute_list_dispatch(compute_list, x_groups, y_groups, z_groups)
					rd.compute_list_end()
					
					pass
				
				elif render_type == RenderType.VELOCITY_SPHERES:
					
					# linearize depth
					var uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(linear_sampler)
					uniform.add_id(fluid_depth_image)
					var fluid_depth_uniform_set := UniformSetCacheRD.get_cache(linearize_depth_shader, 0, [uniform])
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
					uniform.binding = 0
					uniform.add_id(fluid_linear_depth_image)
					var fluid_linear_depth_uniform_set := UniformSetCacheRD.get_cache(linearize_depth_shader, 1, [uniform])
					
					
					var matrices_uniform : RDUniform = RDUniform.new()
					matrices_uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_UNIFORM_BUFFER
					matrices_uniform.binding = 0
					matrices_uniform.add_id(inv_projection_matrix_uniform_buffer)
					var matrices_uniform_set: RID = UniformSetCacheRD.get_cache(linearize_depth_shader, 3, [matrices_uniform])
					
					var push_constant : PackedFloat32Array = PackedFloat32Array()
					push_constant.push_back(render_size.x)
					push_constant.push_back(render_size.y)
					push_constant.push_back(0.0)
					push_constant.push_back(0.0)
					
					@warning_ignore("integer_division")
					var x_groups := (render_size.x - 1) / 8 + 1
					@warning_ignore("integer_division")
					var y_groups := (render_size.y - 1) / 8 + 1
					var z_groups := 1
					
					var compute_list := rd.compute_list_begin()
					rd.compute_list_bind_compute_pipeline(compute_list, linearize_depth_pipeline)
					rd.compute_list_bind_uniform_set(compute_list, fluid_depth_uniform_set, 0)
					rd.compute_list_bind_uniform_set(compute_list, fluid_linear_depth_uniform_set, 1)
					rd.compute_list_bind_uniform_set(compute_list, matrices_uniform_set, 3)
					rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
					rd.compute_list_dispatch(compute_list, x_groups, y_groups, z_groups)
					rd.compute_list_end()
					
					# get sphere color render
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(linear_sampler)
					uniform.add_id(fluid_linear_depth_image)
					fluid_linear_depth_uniform_set = UniformSetCacheRD.get_cache(fluid_render_velocity_spheres_shader, 0, [uniform])
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(linear_sampler)
					uniform.add_id(fluid_color_image)
					var fluid_color_uniform_set = UniformSetCacheRD.get_cache(fluid_render_velocity_spheres_shader, 1, [uniform])
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(linear_sampler)
					uniform.add_id(depth_image)
					var bg_depth_uniform_set := UniformSetCacheRD.get_cache(fluid_render_velocity_spheres_shader, 2, [uniform])
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
					uniform.binding = 0
					uniform.add_id(temp_color_image)
					var temp_color_uniform_set := UniformSetCacheRD.get_cache(fluid_render_velocity_spheres_shader, 3, [uniform])
					
					matrices_uniform = RDUniform.new()
					matrices_uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_UNIFORM_BUFFER
					matrices_uniform.binding = 0
					matrices_uniform.add_id(inv_projection_projection_matrix_uniform_buffer)
					matrices_uniform_set = UniformSetCacheRD.get_cache(fluid_render_velocity_spheres_shader, 4, [matrices_uniform])
					
					compute_list = rd.compute_list_begin()
					rd.compute_list_bind_compute_pipeline(compute_list, fluid_render_velocity_spheres_pipeline)
					rd.compute_list_bind_uniform_set(compute_list, fluid_linear_depth_uniform_set, 0)
					rd.compute_list_bind_uniform_set(compute_list, fluid_color_uniform_set, 1)
					rd.compute_list_bind_uniform_set(compute_list, bg_depth_uniform_set, 2)
					rd.compute_list_bind_uniform_set(compute_list, temp_color_uniform_set, 3)
					rd.compute_list_bind_uniform_set(compute_list, matrices_uniform_set, 4)
					rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
					rd.compute_list_dispatch(compute_list, x_groups, y_groups, z_groups)
					rd.compute_list_end()
					
					# copy color to main 
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(nearest_sampler)
					uniform.add_id(temp_color_image)
					temp_color_uniform_set = UniformSetCacheRD.get_cache(copy_texture_shader, 0, [uniform])
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
					uniform.binding = 0
					uniform.add_id(color_image)
					var color_uniform_set = UniformSetCacheRD.get_cache(copy_texture_shader, 1, [uniform])
					
					push_constant = PackedFloat32Array()
					push_constant.append(render_size.x)
					push_constant.append(render_size.y)
					push_constant.append(0.0) # color_threshold
					push_constant.append(0.0)
					
					compute_list = rd.compute_list_begin()
					rd.compute_list_bind_compute_pipeline(compute_list, copy_texture_pipeline)
					rd.compute_list_bind_uniform_set(compute_list, temp_color_uniform_set, 0)
					rd.compute_list_bind_uniform_set(compute_list, color_uniform_set, 1)
					rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
					rd.compute_list_dispatch(compute_list, x_groups, y_groups, z_groups)
					rd.compute_list_end()
					
					pass
				
				if debug_draw_depth:
				
					# Create uniform sets, these will be cached, the cache will be cleared if our viewports configuration is changed.
					var uniform := RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(linear_sampler)
					uniform.add_id(fluid_depth_image)
					var fluid_depth_uniform_set := UniformSetCacheRD.get_cache(particle_depth_visualiser_shader, 0, [uniform])
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_IMAGE
					uniform.binding = 0
					uniform.add_id(color_image)
					var color_uniform_set := UniformSetCacheRD.get_cache(particle_depth_visualiser_shader, 1, [uniform])
					
					
					uniform = RDUniform.new()
					uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_SAMPLER_WITH_TEXTURE
					uniform.binding = 0
					uniform.add_id(linear_sampler)
					uniform.add_id(fluid_linear_depth_image)
					var fluid_linear_depth_uniform_set := UniformSetCacheRD.get_cache(particle_depth_visualiser_shader, 2, [uniform])
					
					var matrices_uniform : RDUniform = RDUniform.new()
					matrices_uniform.uniform_type = RenderingDevice.UNIFORM_TYPE_UNIFORM_BUFFER
					matrices_uniform.binding = 0
					matrices_uniform.add_id(inv_projection_projection_matrix_uniform_buffer)
					var matrices_uniform_set: RID = UniformSetCacheRD.get_cache(particle_depth_visualiser_shader, 3, [matrices_uniform])
					
					# Create push constant.
					# Must be aligned to 16 bytes and be in the same order as defined in the shader.
					# We don't have structures (yet) so we need to build our push constant "the hard way"...
					var push_constant : PackedFloat32Array = PackedFloat32Array()
					push_constant.push_back(Global.particle_pos_texture_width)
					push_constant.push_back(particle_sphere_radius)
					push_constant.push_back(render_size.x)
					push_constant.push_back(render_size.y)
					push_constant.push_back(depth_divisor)
					push_constant.push_back(0.0)
					push_constant.push_back(0.0)
					push_constant.push_back(0.0)
					
					#@warning_ignore("integer_division")
					#var x_groups : int = (Global.particle_count - 1) / 128 + 1
					#var y_groups : int = 1
					#var z_groups : int = 1
					@warning_ignore("integer_division")
					var x_groups := (render_size.x - 1) / 8 + 1
					@warning_ignore("integer_division")
					var y_groups := (render_size.y - 1) / 8 + 1
					var z_groups := 1
					
					
					# Run our compute shader.
					var compute_list := rd.compute_list_begin()
					rd.compute_list_bind_compute_pipeline(compute_list, particle_depth_visualiser_pipeline)
					rd.compute_list_bind_uniform_set(compute_list, fluid_depth_uniform_set, 0)
					rd.compute_list_bind_uniform_set(compute_list, color_uniform_set, 1)
					rd.compute_list_bind_uniform_set(compute_list, fluid_linear_depth_uniform_set, 2)
					rd.compute_list_bind_uniform_set(compute_list, matrices_uniform_set, 3)
					rd.compute_list_set_push_constant(compute_list, push_constant.to_byte_array(), push_constant.size() * 4)
					rd.compute_list_dispatch(compute_list, x_groups, y_groups, z_groups)
					rd.compute_list_end()
				
#endregion
