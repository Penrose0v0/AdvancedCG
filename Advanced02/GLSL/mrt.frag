#version 130

// TODO: define "in" variables
in vec3 vVertexNormal;
in vec4 cameraCor; 
in float depth; 

out vec4 fragColors[3];

void main()
{	
	// 0v0
	vec3 tmp = 0.5 * vVertexNormal + 0.5; 

	fragColors[0] = cameraCor;
	fragColors[1] = vec4(tmp, 1.0); ;
	fragColors[2] = vec4(depth, depth, depth, 1); 
}
