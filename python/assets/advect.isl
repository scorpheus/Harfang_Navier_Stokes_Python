in {
	float size_pixel;
	float dt0;
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
			vec3 uvdens_in = texture2D(tex_xydens_in, v_uv).xyz;		
			float i = v_uv.x / size_pixel;
			float j = v_uv.y / size_pixel;

			float x = i - dt0 * uvdens_in.x;
			float y = j - dt0 * uvdens_in.y;
			if(x < 0.5f)
				x = 0.5f;
			if(x > N + 0.5f)
			   x = N + 0.5f;
			float i0 = floor(x);
			float i1 = i0 + 1.0f;
			if(y < 0.5)
				y = 0.5;
			if(y > N + 0.5f)
			   y = N + 0.5f;
			float j0 = floor(y);
			float j1 = j0 + 1.0f;
			float s1 = x - i0;
			float s0 = 1.0f - s1;
			float t1 = y - j0;
			float t0 = 1.0f - t1;
			%out.color% = vec4(s0 * (t0 * GetNeighborPixelColor(v_uv, i0 - i, j0 - j) + t1 * GetNeighborPixelColor(v_uv, i0 - i, j1 - j)) + s1 * (t0 * GetNeighborPixelColor(v_uv, i1 - i, j0 - j) + t1 * GetNeighborPixelColor(v_uv, i1 - i, j1 - j)), 1);
        %}
    }
}
