in {
    tex2D tex_xydens [wrap:clamp];
	vec2 size_pixel;
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
            vec3 GetNeighborPixelColor(vec2 v_uv, int x, int y){return texture2D(tex_xydens, (v_uv + size_pixel * vec2(x, y))).xyz;}
		%}

        source %{
			
           // vec4 t = vec4(GetNeighborPixelColor(v_uv, -1, -1), 1);
            %out.color% = vec4(10, 0, 0, 1);
        %}
    }
}
