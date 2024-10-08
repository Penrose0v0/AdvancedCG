#version 130

in vec3 ePosition;
in vec3 eNormal;

out vec4 fragColor;

uniform float roughness;
uniform vec3 fresnel0;
uniform vec3 eLightDir;
uniform vec3 lightColor;
uniform vec3 diffuseCoeff;
uniform vec3 ambient;

//#define PI 3.14159265358979323846264338327950288

void main()
{
	vec3 eNormalNormalized = normalize(eNormal);
	float dotDiffuse = max(dot(eNormalNormalized, eLightDir), 0.0);
	vec3 diffuseColor = dotDiffuse * lightColor * diffuseCoeff;

	vec3 specularColor = vec3(0, 0, 0);
	if (dotDiffuse > 0.0)
	{
		vec3 viewingDir = normalize(-ePosition);
		vec3 halfVec = normalize(eLightDir + viewingDir);

		float NH = dot(eNormalNormalized, halfVec);
		float VH = dot(viewingDir, halfVec);
		float NV = dot(eNormalNormalized, viewingDir);
		float NL = dot(eNormalNormalized, eLightDir);

		float m2 = roughness * roughness;
		float D = (1.0/(4*m2*NH*NH*NH*NH))*exp((NH*NH-1)/(m2*NH*NH));

		float G = min(1, min((2*NH*NV)/VH, (2*NH*NL)/VH));

		float fresnelCoeff = 1.0 - NV;
		float fresnelCoeff2 = fresnelCoeff * fresnelCoeff;
		vec3 F = fresnel0 + (vec3(1.0) - fresnel0) * (fresnelCoeff2 * fresnelCoeff2 * fresnelCoeff);

		specularColor = D * G * F / (4.0 * NV) * lightColor;
	}

	fragColor = vec4(specularColor + diffuseColor + ambient, 1);
}
