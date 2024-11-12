#version 130

in vec3 eNormal;

out vec4 fragColor;

uniform vec3 eLightDir;
uniform vec3 lightColor;
uniform vec3 diffuseCoeff;
uniform vec3 ambient;

void main()
{
	float dotProd = max(dot(normalize(eNormal), eLightDir), 0.0);
	fragColor = vec4(dotProd * lightColor * diffuseCoeff + ambient, 1);
}
