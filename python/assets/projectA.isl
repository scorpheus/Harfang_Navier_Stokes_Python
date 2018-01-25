in {
	float size_pixel;
	float N;
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
		
		global %{
            vec3 GetNeighborPixelColor(vec2 v_uv, float x, float y){return texture2D(tex_xydens_in, (v_uv + size_pixel * vec2(x, y))).xyz;}
		%}

        source %{
			
			vec3 uvdens_out = texture2D(tex_xydens_out, v_uv).xyz;
			float u = GetNeighborPixelColor(v_uv, 1, 0).x - GetNeighborPixelColor(v_uv, -1, 0).x; 
			float v = GetNeighborPixelColor(v_uv, 1, 0).y - GetNeighborPixelColor(v_uv, -1, 0).y; 

			float y = -0.5 * (u + v) / N;
			%out.color% = vec4(0, y, uvdens_out.z, 1);
        %}
    }
}
