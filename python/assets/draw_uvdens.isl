in {
	float dt;
    tex2D tex_xydens_in [wrap:clamp];
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
			%out.color% = vec4(texture2D(tex_xydens_in, v_uv).xyz/1.0f, 1);
        %}
    }
}
