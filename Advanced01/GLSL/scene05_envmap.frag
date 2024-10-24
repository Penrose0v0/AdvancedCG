﻿#version 130

#define PI 3.141592653589793

in vec3 vWorldEyeDir;
in vec3 vWorldNormal;

out vec4 fragColor;

// TODO: uncomment these lines
uniform sampler2D envmap;

float atan2(in float y, in float x)
{
    return x == 0.0 ? sign(y)*PI/2 : atan(y, x);
}

void main()
{
	// TODO: write an appropriate code here
	vec3 reflectedDir = reflect(vWorldEyeDir, vWorldNormal);

	float theta = atan2(reflectedDir.y, reflectedDir.x); 
	// float phi = reflectedDir.x > 0 ? asin(reflectedDir.z) : asin(reflectedDir.z) + sign(reflectedDir.z) * PI / 2; 
	float phi = asin(reflectedDir.z); 

	float latitude = (theta + PI) / (2.0 * PI);
    float longitude = (phi + PI / 2.0) / PI; 

	fragColor = texture2D(envmap, vec2(longitude, latitude));
}
