@tool
class_name MainCamera
extends Camera3D
## This is a custom node meant to hold several child [Camera3D]s and their connected [SubViewport]s.
## It ensures that the properties of the cameras are all synced up with the properties of this main
## camera. It also ensures that the sizes of the subviewports all match the viewport which is a
## parent of this camera. This means that the textures from the subviewports will be the same
## resolution as the main viewport.



# PUBLIC PROPERTIES

## Holds an array of properties that are unique to [Camera3D]s.
static var camera_3d_unique_properties: Array[StringName] = get_camera_3d_unique_properties()



# PRIVATE PROPERTIES

# Camera3D properties which should not be copied to the child cameras.
var _excluded_properties: Array[StringName] = [&"cull_mask", &"current", &"environment", &"compositor"]

var ssfr: ScreenSpaceFluidRendering

# Keep references to various child nodes.
@onready var _viewports: Array[Node] = find_children("*", "SubViewport")
@onready var _cameras: Array[Node] = find_children("*", "Camera3D")

# Keep a reference to the viewport which is above this camera.
@onready var _parent_viewport: Viewport = get_viewport()

# PUBLIC METHODS

## Returns a list of property names unique to [Camera3D]s.
static func get_camera_3d_unique_properties() -> Array[StringName]:
	var child_camera_settable_properties: Array[StringName] = []
	for prop: Dictionary in ClassDB.class_get_property_list("Camera3D", true):
		child_camera_settable_properties.append(prop["name"])
	return child_camera_settable_properties


func set_ssfr_render_type(new_render_type: ScreenSpaceFluidRendering.RenderType) -> void:
	ssfr.render_type = new_render_type

# PRIVATE METHODS

# Called when the node enters the scene tree for the first time.
func _ready() -> void:
	
	if not Engine.is_editor_hint():
		Global.current_camera = self
	
	# find ssfr effect and give it the first directional light's facing direction
	for effect in self.compositor.compositor_effects:
		if effect is ScreenSpaceFluidRendering:
			ssfr = effect
	
	if not Engine.is_editor_hint():
		var directional_lights_array := self.get_parent().find_children("*", "DirectionalLight3D")
		ssfr.light_world_dir = -directional_lights_array[0].global_transform.basis.z
		print("directional light basis -z: ", -directional_lights_array[0].global_transform.basis.z)
		pass
	
	# Connect to changes in the parent viewport
	if is_instance_valid(_parent_viewport):
		_parent_viewport.size_changed.connect(_on_parent_viewport_size_changed)
		_on_parent_viewport_size_changed.call_deferred()

# Called whenever the parent viewport's size changes.
func _on_parent_viewport_size_changed() -> void:
	if !Engine.is_editor_hint() and is_instance_valid(_parent_viewport):
		# Update the size of any child viewport nodes
		for viewport: Node in _viewports:
			if viewport is SubViewport:
				viewport.size = _parent_viewport.size

# Overriding what happens when a property is set.
# Keeps the depth and normals camera in sync with the main camera.
func _set(property: StringName, value: Variant) -> bool:
	# This branch sets both the main camera and the child cameras
	if property in camera_3d_unique_properties and property not in _excluded_properties:
		for camera: Node in _cameras:
			if camera is Camera3D:
				camera.set(property, value)
		return false
	# This branch only sets the main camera
	else:
		return false


## Flying camera part

const SPEED = 200.0

@export var flying_enabled : bool = false

var current_rotation : Vector2 = Vector2.ZERO 

func _unhandled_input(event: InputEvent) -> void:
	
	if (event is InputEventMouseButton):
		if (event.button_index == MOUSE_BUTTON_RIGHT):
			if (event.pressed == true):
				flying_enabled = true
			else:
				flying_enabled = false
	
	# mouse look
	if (event is InputEventMouseMotion) and flying_enabled:
		current_rotation -= event.relative * 0.2
		
		if (abs(current_rotation.x) > 360.0):
			current_rotation.x = 0.0;
		
		if (abs(current_rotation.y) > 89.9):
			current_rotation.y = sign(current_rotation.y) * 89.9;
		
		pass

func _process(delta: float) -> void:
	
	if flying_enabled:
		# mouse look
		self.global_rotation_degrees = Vector3(current_rotation.y, current_rotation.x, 0.0)
		
		# move
		var input_dir : Vector2 = Input.get_vector("move_left", "move_right", "move_forward", "move_backward")
		var direction : Vector3 = (self.transform.basis * Vector3(input_dir.x, 0, input_dir.y)).normalized()
		
		self.position += direction * SPEED * delta
	
