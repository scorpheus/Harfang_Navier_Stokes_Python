in {
	float dt;
    tex2D tex_xydens_in [wrap:clamp];
    tex2D tex_xydens_out [wrap:clamp];
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
			vec3 uvdens_in = texture2D(tex_xydens_in, v_uv).xyz;
			vec3 uvdens_out = texture2D(tex_xydens_out, v_uv).xyz;
			%out.color% = vec4(uvdens_out.xyz + dt * uvdens_in.xyz, 1);
        %}
    }
}
