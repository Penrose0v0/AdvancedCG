#version 130

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
	float theta = atan2(vWorldNormal.y, vWorldNormal.x); 
	float sai = atan2(vWorldNormal.z, vWorldNormal.x); 
	fragColor = texture2D(envmap, vec2(sai, theta));
}
