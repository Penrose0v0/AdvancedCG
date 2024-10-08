#version 150 core

in vec3 ePosition;
in vec3 eNormal;
in vec3 vColor;

out vec4 fragColor;

void main()
{
	float dotProd = abs(dot(normalize(ePosition), normalize(eNormal)));
	fragColor = vec4(dotProd * vColor, 1);
}
