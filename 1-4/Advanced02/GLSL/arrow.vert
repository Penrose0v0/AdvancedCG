#version 130

in vec4 vertexPosition;
in vec3 vertexNormal;
in vec3 vertexColor;
in vec2 sphericalCoord;

out vec3 ePosition;
out vec3 eNormal;
out vec3 vColor;

uniform float scale;
uniform float radius;
//uniform vec2 sphericalCoord;
uniform mat4 modelViewMatrix;
uniform mat4 projMatrix;

mat3 sphericalCoord2rotation(in vec2 sphericalCoord)
{
	float cosPhi = cos(sphericalCoord.x), sinPhi = sin(sphericalCoord.x);
	float cosTheta = cos(sphericalCoord.y), sinTheta = sin(sphericalCoord.y);
	return mat3(cosPhi, 0.0, -sinPhi, sinPhi * sinTheta, cosTheta, cosPhi * sinTheta, sinPhi * cosTheta, -sinTheta, cosPhi * cosTheta);
}

void main()
{
	mat3 R = sphericalCoord2rotation(sphericalCoord);
	vec3 modelTransformed = R * (scale * vec3(1,-1,1) * vertexPosition.xyz + vec3(0, radius, 0));
	ePosition = (modelViewMatrix * vec4(modelTransformed, 1)).xyz;
	eNormal = mat3(modelViewMatrix) * R * vertexNormal;
	vColor = vertexColor;
	gl_Position = projMatrix * vec4(ePosition, 1.0);
}