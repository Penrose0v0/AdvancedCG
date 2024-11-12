#version 150 core

#define PI 3.141592653589793

in vec3 vWorldEyeDir;
in vec3 vWorldNormal;

out vec4 fragColor;

// TODO: uncomment these lines
//uniform sampler2D envmap;

float atan2(in float y, in float x)
{
    return x == 0.0 ? sign(y)*PI/2 : atan(y, x);
}

void main()
{
	// TODO: write an appropriate code here
	fragColor = vec4(0, 1, 0, 1);
}
