in {
	tex2D array_tex [wrap:clamp];
	float width;
	float height;
	float depth;
	float size_array;
}

variant {
	vertex {
		source %{
			%out.position% = vec4(vPosition, 1.0);
		%}
	}

	pixel {
		global %{
			bool hitbox(vec3 origin, vec3 direction, vec3 m1, vec3 m2, out float tmin, out float tmax) 
			 {
			   float tymin, tymax, tzmin, tzmax; 
			   float flag = 1.0; 
 
				if (direction.x >= 0) 
				{
				   tmin = (m1.x - origin.x) / direction.x;
					 tmax = (m2.x - origin.x) / direction.x;
				}
				else 
				{
				   tmin = (m2.x - origin.x) / direction.x;
				   tmax = (m1.x - origin.x) / direction.x;
				}
				if (direction.y >= 0) 
				{
				   tymin = (m1.y - origin.y) / direction.y; 
				   tymax = (m2.y - origin.y) / direction.y; 
				}
				else 
				{
				   tymin = (m2.y - origin.y) / direction.y; 
				   tymax = (m1.y - origin.y) / direction.y; 
				}
     
				if ((tmin > tymax) || (tymin > tmax)) flag = -1.0; 
				if (tymin > tmin) tmin = tymin; 
				if (tymax < tmax) tmax = tymax; 
      
				if (direction.z >= 0) 
				{
				   tzmin = (m1.z - origin.z) / direction.z; 
				   tzmax = (m2.z - origin.z) / direction.z; 
				}
				else 
				{
				   tzmin = (m2.z - origin.z) / direction.z; 
				   tzmax = (m1.z - origin.z) / direction.z; 
				}
				if ((tmin > tzmax) || (tzmin > tmax)) flag = -1.0; 
				if (tzmin > tmin) tmin = tzmin; 
				if (tzmax < tmax) tmax = tzmax; 
      
				return (flag > 0); 
			 }

			void intersection_distances_no_if(vec3 origin, vec3 direction, vec3 aabb[2], out float tmin, out float tmax){
				vec3 inv_direction = vec3(1.0) / direction;
				int sign_0 = inv_direction.x < 0 ? 1 : 0;
				int sign_1 = inv_direction.y < 0 ? 1 : 0;
				int sign_2 = inv_direction.z < 0 ? 1 : 0;

				float tymin, tymax, tzmin, tzmax;
				tmin = (aabb[sign_0].x - origin.x) * inv_direction.x;
				tmax = (aabb[1-sign_0].x - origin.x) * inv_direction.x;
				tymin = (aabb[sign_1].y - origin.y) * inv_direction.y;
				tymax = (aabb[1-sign_1].y - origin.y) * inv_direction.y;
				tzmin = (aabb[sign_2].z - origin.z) * inv_direction.z;
				tzmax = (aabb[1-sign_2].z - origin.z) * inv_direction.z;
				tmin = max(max(tmin, tymin), tzmin);
				tmax = min(min(tmax, tymax), tzmax);
			}
			
			const int NUM_STEPS = 100;
			float heightMapTracing(vec3 ori, vec3 dir, float max_dist, out vec3 p) {
				float tot = 0.f;
				float substep = 0.01f;
				//float max_num_step = min(max_dist / substep, NUM_STEPS);
				float max_num_step = max_dist / substep;
				float x, y, z;
				for(int i = 0; i < max_num_step; i++) {
					p = ori + dir * i * substep;
					x = floor((p.x + 0.5f)*width);
					y = floor((p.y + 0.5f)*width);
					z = floor((p.z + 0.5f)*width);

					float u = (x + width * y + width * height * z) / size_array;

					float value = texture2D(array_tex, vec2(u, 0.5f)).x;
					tot += value;
				}
				return tot;
			}
		%}
		source %{
			vec2 uv = %in.fragcoord%.xy * vInverseInternalResolution.xy;
			uv = uv * vec2(2.0) - vec2(1.0);

			// ratio
			uv.y *= (vViewport.w / vViewport.z);			

			// ray
			vec3 ori = vViewPosition.xyz;
			vec3 dir = normalize(vec3(uv.xy, 2.5));
			dir = normalize((vec4(dir, 1.0) * vViewMatrix).xyz);
			
			float tmin = 0, tmax = 10000;
			vec3 aabb[2];
			vec3 scale = vec3(1.0, height/width, depth/width);
			aabb[0] = vec3(-0.5f, -0.5f, -0.5f);
			aabb[1] = scale - vec3(0.5f, 0.5f, 0.5f);		
			if(!hitbox(ori, dir, aabb[0], aabb[1], tmin, tmax))
				discard;

			//intersection_distances_no_if(ori, dir, aabb, tmin, tmax);

			// tracing
			vec3 p;
	//		tmin = max(tmin, 0.0);
	//		tmax = max(tmax, 0.0);

		//	if(tmax <= tmin)
		//		discard;
			
			float tot = heightMapTracing(ori + dir * tmin, dir, tmax - tmin, p);
			if(tot > 0.0)
				%out.color% = vec4(tot, 0.0, tot, 1.0);
			else
				%out.color% = vec4(0, 0.5, 0.0, 0.0);
				
		%}
	}
}
