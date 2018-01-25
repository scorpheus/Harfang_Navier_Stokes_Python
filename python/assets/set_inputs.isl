in {
	vec2 size_pixel;
	vec3 left_button_pressed;
}

variant {
    vertex {
        out { vec2 v_uv; }

        source %{
            v_uv = vUV0;
			%out.position% = vec4(vPosition, 1.0);
        %}
    }

    pixel {
        in { vec2 v_uv; }
		
        source %{
			vec4 xydens = vec4(0, 0, 0, 1);

			// get the pressed mouse bouton and put density on it
			if(left_button_pressed.z == 1.0 && abs(left_button_pressed.x - v_uv.x) <= size_pixel.x*5 && abs(left_button_pressed.y - v_uv.y) <= size_pixel.y*5 )
				xydens.z = 100.0f;
			// add force on some	
			if(v_uv.x <= 0.1 && v_uv.y > 0.9 )
				xydens.xy = vec2(5.0f, 5.0f);

            %out.color% = xydens;
        %}
    }
}
