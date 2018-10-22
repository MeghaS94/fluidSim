#version 330 core
out vec4 vertexColor;
layout (location = 0) in vec3 aPos; 
layout (location = 1) in vec3 densitySource; /* source of density */
uniform mat4 transform;
uniform mat4 projection;

void main()
{
	gl_Position = projection*transform * vec4(aPos, 1.0);
    // color should actually be densitySource.x, densitySource.y, densitySource.z
    vertexColor = vec4(0.5, 0.5, 0.5, 1.0f);
}

