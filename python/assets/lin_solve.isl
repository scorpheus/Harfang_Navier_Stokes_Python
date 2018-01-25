in {
	float size_pixel;
	float a;
	float c;
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
            vec3 GetNeighborPixelColor(vec2 v_uv, int x, int y){return texture2D(tex_xydens_out, (v_uv + size_pixel * vec2(x, y))).xyz;}
		%}

        source %{
			vec3 uvdens_in = texture2D(tex_xydens_in, v_uv).xyz;			
		
			%out.color% = vec4((uvdens_in + a * (GetNeighborPixelColor(v_uv, -1, 0) + GetNeighborPixelColor(v_uv, 1, 0) + GetNeighborPixelColor(v_uv, 0, -1) + GetNeighborPixelColor(v_uv, 0, 1)))/ c, 1);
        %}
    }
}
