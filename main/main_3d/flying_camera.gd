extends Camera3D

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
		current_rotation -= event.relative
		
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
	
