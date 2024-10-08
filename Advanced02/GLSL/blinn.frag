#version 130

in vec3 ePosition;
in vec3 eNormal;

out vec4 fragColor;

void main()
{
	float dotProd = abs(dot(normalize(ePosition), normalize(eNormal)));
	fragColor = vec4(vec3(dotProd), 1);
}
