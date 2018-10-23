/*
Date : 22/08/18
Command to compile : g++ -std=c++11 main.cc Shader.cc glad.c -lglfw -lGL -lX11 -lpthread -lXrandr -lXi -ldl -o hello
g++ -std=c++11 src/main.cc src/Shader.cc src/glad.c src/solverUtils.cc -lglfw3 -lGLU -lGL -lX11 -lXxf86vm -lpthread -lXrandr -lXi -ldl -lXinerama -lXcursor -lrt -lm -o hello
*/

#include <glad/glad.h>	// glad manages hardware specific function pointers for openGL, gives address of extensions etc
#include <GLFW/glfw3.h>
#include <iostream>
#include <vector>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>	

#include "../includes/Shader.h"	
#include "../includes/solverUtils.h"

using namespace std;

//macros 
//array value accessor
#define IX(i,j) ((i)+(N+2)*(j))
//resolution
#define N 100

// grid of density sources
float *dens_prev;
float *dens;
float *u, *v, *u_prev, *v_prev;

int WINDOW_SIZE_WIDTH = 600;
int WINDOW_SIZE_HEIGHT = 600;

// vertices for a grid
float gridPoints[4*(N+2)*(N+2)];

// grid size
#define Size (N+2)*(N+2)

float vel_pos[6*(N+2)*(N+2)];
float dens_grid[4*Size];
float dens_quads[36*N*N];

static float dt = 0.004f;
static float source = 20.0f;
// rate of diffusion
static float diff = 0.0004f;

// callback function when window is resized to adjust viewport
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    glViewport(0, 0, width, height);
	WINDOW_SIZE_WIDTH = width;
	WINDOW_SIZE_HEIGHT = height;
}  

void processInput(GLFWwindow *window)
{
    if(glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
	{
        glfwSetWindowShouldClose(window, true);
	}
	if (glfwGetKey(window, GLFW_KEY_E) == GLFW_PRESS)
	{
		std::cout << "key e pressed" << std::endl;
	}
	double xpos, ypos;
	glfwGetCursorPos(window, &xpos, &ypos);
}

// remap values from mouse coordinates to screen/ortho projection coordinates
float valueRemapping(float value, float low1, float high1, float low2, float high2)
{
	return low2 + (value - low1) * (high2 - low2) / (high1 - low1);
}	

void mouse_button_callback(GLFWwindow* window, int button, int action, int mods)
{
    if (button == GLFW_MOUSE_BUTTON_RIGHT && action == GLFW_PRESS)
	{
		cout << "right mouse button" << endl;
		double xpos, ypos;
		glfwGetCursorPos(window, &xpos, &ypos);
		glm::vec2 result;
		result.x = valueRemapping(xpos, 0, WINDOW_SIZE_WIDTH, 0, N);
		result.y = valueRemapping(ypos, 0, WINDOW_SIZE_HEIGHT, 0, N);
		// std::cout << xpos/double(WINDOW_SIZE) << ", " << 1.0 -ypos/double(WINDOW_SIZE)  << std::endl;
		std::cout << result.x << ", " << result.y << endl;
		dens_prev[IX(int(result.x), int((N+2) - result.y))] = 20000.0;
		
		int max = 10.0 , min = -10.0;
		// u_prev[IX(int(result.x), int((N+2) - result.y))] = (float(rand()%(max-min + 1) + min)/10.0f ) * 200.0f;
		// v_prev[IX(int(result.x), int((N+2) - result.y))] = (float(rand()%(max-min + 1) + min)/10.0f) * 200.0f;
		
		std::cout << "density sources value at " << result.x << ", "<< (N+2) - result.y << " : " << dens_prev[IX(int(result.x), int((N+2) - result.y))] << endl;
		// int index = 0;
		// for (int i=0;i<N+2;i++)
		// {
		// 	for (int j=0;j<N+2;j++)
		// 	{
		// 		gridPoints[index] = float(i);///float(N);
		// 		gridPoints[index+1] = float(j);///float(N);
		// 		gridPoints[index+2] = 0.0;  
		// 		gridPoints[index+3] += dens_prev[IX(i,j)]*20;
		// 		// std::cout << gridPoints[index]<< std::endl;
		// 		index = index + 4;
		// 	}
		// }
	}
	if (button == GLFW_MOUSE_BUTTON_LEFT && action == GLFW_PRESS)
	{
		cout << "left mouse button" << endl;
	}
}

struct vec3f{
	double x;
	double y;
	double z;
};

std::vector<vec3f> cursor_pos;


void addVelocitySource(std::vector<vec3f> cursor_pos, int mode)
{
	// if mode == 1, then add mouse direction , else add random
	// modify u_prev and v_prev with the cursor movement
	if (mode == 2)
	{
		int numElements = cursor_pos.size();
		if (numElements % 2 != 0)
		{
			// the num elements are not even - pairs of 2
			cursor_pos.pop_back();
		}

			for (int i=0;i<cursor_pos.size();i = i+2)
			{
				double x1 = cursor_pos[i].x;
				double y1 = cursor_pos[i].y;

				double x2 = cursor_pos[i+1].x;
				double y2 = cursor_pos[i+1].y;

				int gridX = int(x1), gridY = int(y1);
				if (abs(int(x1) - x1) > 0.5)
				{
					gridX = int(x1) + 1;
				}
				if (abs(int(y1) - y1) > 0.5)
				{
					gridY = int(y1) + 1;
				}
				
				double dx = x2 - x1;
				double dy = y2 - y1;
				u_prev[IX(gridX, gridY)] = dx*10.0f;
				v_prev[IX(gridX, gridY)] = dy*10.0f; 
				// std::cout << dx << ", " << dy << std::endl;
			}
	}
}

static void cursor_position_callback(GLFWwindow* window, double xpos, double ypos)
{
	if (glfwGetMouseButton(window, GLFW_MOUSE_BUTTON_LEFT) == GLFW_RELEASE) 
  	{
		return;
  	}
	
	int index_x = int(valueRemapping(xpos, 0, WINDOW_SIZE_WIDTH, 0, N));
	int index_y = int((N+2) - valueRemapping(ypos, 0, WINDOW_SIZE_WIDTH, 0, N));
	double x = valueRemapping(xpos, 0, WINDOW_SIZE_WIDTH, 0, N);
	double y = (N+2) - valueRemapping(ypos, 0, WINDOW_SIZE_WIDTH, 0, N);
	cursor_pos.push_back(vec3f());
	cursor_pos[cursor_pos.size()-1].x  = x;
	cursor_pos[cursor_pos.size()-1].y  = y;
	cursor_pos[cursor_pos.size()-1].z  = 0.0f; 
	
	addVelocitySource(cursor_pos, 2);
}
// not used yet
void drawGrid(float* gridValues, unsigned int VBO, unsigned int VAO, int shaderID)
{
	glBufferData(GL_ARRAY_BUFFER, sizeof(gridValues), gridValues, GL_STREAM_DRAW);
	
	GLint posAttrib = glGetAttribLocation(shaderID, "aPos");
	std::cout << posAttrib << std::endl;
	// iterpreting data from buffer 
	glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE, 0, (void*)0);
	// specify index of the vertex attribute to be enabled. --> 0
    glEnableVertexAttribArray(0);

	GLint densAttrib = glGetAttribLocation(shaderID, "densitySource");
	std::cout << densAttrib << std::endl;
	// iterpreting data from buffer 
	glVertexAttribPointer(densAttrib, 3, GL_FLOAT, GL_FALSE, 4 * sizeof(float), (void*)(3* sizeof(float)));
	// specify index of the vertex attribute to be enabled. --> 0
    glEnableVertexAttribArray(1);
}
 
void draw_density()
{
    float x;
    float y;
    float d00;
    float d01;
    float d10;
    float d11;

    int rowSize = N+2;
    int colSize = N+2;

	int index = 0;
	

        for(int i=1; i<=N; i++)
        {
            x = (float)i;
            for(int j=1; j<=N; j++)
            {
                y = (float)j;

                d00 = dens_prev[IX(i, j)];
                d01 = dens_prev[IX(i, j+1)];
                d10 = dens_prev[IX(i+1, j)];
                d11 = dens_prev[IX(i+1, j+1)];
				// std::cout << d00 << ", " << d01 << ", " << d11 << ", " << d10 << std::endl;
				// std::cout << "x, y | x, y+1 | x+1, y | x+1 , y+1 : " << 
				// 			x<< "," << y << " | " << x << ","<< y+1 << " | " << x+1 << "," << y << " | " <<
																							// x+1 << "," << y+1  << std::endl;
				// v0
				dens_quads[index] =   x;
				dens_quads[index+1] =  y;
				dens_quads[index+2] =  0.0f;

				dens_quads[index+3] =  1.0f - d00;
				dens_quads[index+4] =   1.0f;
				dens_quads[index+5] =   1.0f - d00;

				// v1
				dens_quads[index+6] =   x +1;
				dens_quads[index+7] =   y ;
				dens_quads[index+8] =   0.0f;

				dens_quads[index+9] =  1.0f - d10;
				dens_quads[index+10] =  1.0f;
				dens_quads[index+11] =  1.0f - d10;
				// v2
				dens_quads[index+12] =   x ;
				dens_quads[index+13] =   y + 1.0f;
				dens_quads[index+14] =   0.0f;

				dens_quads[index+15] =  1.0f - d01;
				dens_quads[index+16] =  1.0f;
				dens_quads[index+17] =  1.0f - d01;

				// v2
				dens_quads[index+18] =   x + 1.0f;
				dens_quads[index+19] =   y ;
				dens_quads[index+20] =   0.0f;
				
				dens_quads[index+21] =  1.0f - d10;
				dens_quads[index+22] =  1.0f;
				dens_quads[index+23] =  1.0f - d10;

				// v3
				dens_quads[index+24] =   x + 1.0f;
				dens_quads[index+25] =  y + 1.0f;
				dens_quads[index+26] =  0.0f;	

				dens_quads[index+27] =  1.0f - d11;
				dens_quads[index+28] =  1.0f;
				dens_quads[index+29] =  1.0f - d11;


				// v0
				dens_quads[index+30] =   x;
				dens_quads[index+31] =  y + 1.0f;
				dens_quads[index+32] =  0.0f;	

				dens_quads[index+33] =  1.0f - d01;
				dens_quads[index+34] =  1.0f;
				dens_quads[index+35] =  1.0f - d01;	

				index += 36;
            }
        }
}

// should this function is static
void allocate_data()
{
	// int size = (M+2)*(N+2)*(O+2);

	u			= (float *) malloc ( Size*sizeof(float) );
	v			= (float *) malloc ( Size*sizeof(float) );
	// w			= (float *) malloc ( size*sizeof(float) );
	u_prev		= (float *) malloc ( Size*sizeof(float) );
	v_prev		= (float *) malloc ( Size*sizeof(float) );
	// w_prev		= (float *) malloc ( size*sizeof(float) );
	dens		= (float *) malloc ( Size*sizeof(float) );	
	dens_prev	= (float *) malloc ( Size*sizeof(float) );


	if ( !dens || !dens_prev ) {
		fprintf ( stderr, "cannot allocate data\n" );
	}
	// initialise grid!
	for (int i=0;i<Size;i++)
	{
		dens[i] = dens_prev[i] = u[i] = u_prev[i] = v[i] = v_prev[i] = 0.0f;
	}
}

// should this function be static ?
void free_data()
{
	for (int i=0 ; i<Size ; i++ ) 
	{
		dens[i] = dens_prev[i] = 0.0f;
	}
}

int main()
{
	glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
	glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);


	GLFWwindow* window = glfwCreateWindow(WINDOW_SIZE_WIDTH, WINDOW_SIZE_HEIGHT, "Fluid Sim", NULL, NULL);
	if (window == NULL)
	{
		std::cout << "Failed to create GLFW window" << std::endl;
		glfwTerminate();
		return -1;
	}
	// tell GLFW to make the context of our window the main context on the current thread.
	glfwMakeContextCurrent(window);

	// window resize
	glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);  

	// initialize GLAD
	if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
	{
    	std::cout << "Failed to initialize GLAD" << std::endl;
    	return -1;
	} 

	 // build and compile our shader program
    // ------------------------------------
    Shader ourShader("vertexShader.vs", "fragmentShader.fs"); 
	 
	// set up vertex data (and buffer(s)) and configure vertex attributes
    // ------------------------------------------------------------------

	allocate_data();

	int index = 0;
	for (int i=0;i<N+2;i++)
	{
		for (int j=0;j<N+2;j++)
		{
			gridPoints[index] = float(i);
			gridPoints[index+1] = float(j);
			gridPoints[index+2] = 0.0;  
			gridPoints[index+3] = dens_prev[IX(i,j)];
			index = index + 4;
		}
	}
	
	// vel grid
	for (int i=0;i<Size;i++)
	{
		u_prev[i] = v_prev[i] = 0.5f; u[i] = v[i] = 0.5f;
	}


	int index1 = 0;
	for (int i=0;i<N+2;i++)
	{
		for (int j=0;j<N+2;j++)
		{
			vel_pos[index1] = float(i);
			vel_pos[index1+1] = float(j);
			vel_pos[index1+2] = 0.0f;
			vel_pos[index1+3] = float(i) + u_prev[IX(i,j)];
			vel_pos[index1+4] = float(j) + v_prev[IX(i,j)];
			vel_pos[index1+5] = 0.0f;
			index1 = index1 + 6;
		}
	}

	int index2 =0;
	for (int i=0;i<N+2;i++)
	{
		for (int j=0;j<N+2;j++)
		{
			dens_grid[index2] = float(i);
			dens_grid[index2+1] = float(j);
			dens_grid[index2+2] = 0.0f;
			dens_grid[index2+3] = dens_prev[IX(i,j)];
			index2 += 4;
		}
	}

	// generate buffers- vertex array object and vertex buffer object
	GLuint vaos[2];
    // unsigned int VBO, VAO, density_VBO;
	glGenVertexArrays(2, vaos);
	GLuint vbos[2];
    glGenBuffers(2, vbos);

	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	// render loop
	while(!glfwWindowShouldClose(window))
	{
		draw_density();

		glBindVertexArray(vaos[0]);
		glBindBuffer(GL_ARRAY_BUFFER, vbos[0]);
		glBufferData(GL_ARRAY_BUFFER, sizeof(vel_pos), vel_pos, GL_STREAM_DRAW);

		GLint velAttrib = glGetAttribLocation(ourShader.ID, "aPos");
		// std::cout << posAttrib << std::endl;
		// iterpreting data from buffer 
		glVertexAttribPointer(velAttrib, 3, GL_FLOAT, GL_FALSE,  3 *sizeof(float), (void*)0);
		// specify index of the vertex attribute to be enabled. --> 0
		glEnableVertexAttribArray(0);

		glBindVertexArray(vaos[1]);
		glBindBuffer(GL_ARRAY_BUFFER, vbos[1]);
		glBufferData(GL_ARRAY_BUFFER, sizeof(dens_quads), dens_quads, GL_STREAM_DRAW);

		GLint densAttrib = glGetAttribLocation(ourShader.ID, "densitySource");
		// iterpreting data from buffer 
		glVertexAttribPointer(densAttrib, 3, GL_FLOAT, GL_FALSE,  6 *sizeof(float), (void*)(3* sizeof(float)));
		// specify index of the vertex attribute to be enabled. --> 0
		glEnableVertexAttribArray(0);

		GLint posAttrib = glGetAttribLocation(ourShader.ID, "aPos");
		// std::cout << posAttrib << std::endl;
		// iterpreting data from buffer 
		glVertexAttribPointer(posAttrib, 3, GL_FLOAT, GL_FALSE,  6 *sizeof(float), (void*)0);
		// specify index of the vertex attribute to be enabled. --> 0
		glEnableVertexAttribArray(1);

		// begin simulation
		// input
		processInput(window);
		glfwSetMouseButtonCallback(window, mouse_button_callback);
		glfwSetCursorPosCallback(window, cursor_position_callback);
		
		// rendering commands
		// ...
		// glClearColor is a state setting function, glClear is a state using function
		glClearColor(0.2f, 0.2f, 0.2f, 1.0f); 
		glClear(GL_COLOR_BUFFER_BIT);

	
		// this shader program/shaders are used for all subsequent rendering calls.
		ourShader.use();

		glm::mat4 trans = glm::mat4(1.0f);
		trans = glm::translate(trans, glm::vec3(-0.5f, -0.5f, 0.0f));
		unsigned int transformMatrixLocation = glGetUniformLocation(ourShader.ID, "transform");
		// arg 2- how many matrices, arg3-> should matrix be transposed
		glUniformMatrix4fv(transformMatrixLocation, 1, GL_FALSE, glm::value_ptr(trans));

		glm::mat4 projection = glm::ortho(-10.0f, 110.0f, -1.0f, 110.0f, -1.0f, 100.0f);

		unsigned int projectionMatrixLocation = glGetUniformLocation(ourShader.ID, "projection");
		// arg 2- how many matrices, arg3-> should matrix be transposed
		glUniformMatrix4fv(projectionMatrixLocation, 1, GL_FALSE, glm::value_ptr(projection));

		glBindVertexArray(vaos[0]);
		// glPointSize(4.0f);
		glLineWidth(1.0f);
	    glDrawArrays(GL_LINES,  0, 2*Size);
		
		glBindVertexArray(vaos[1]); 
	    glDrawArrays(GL_TRIANGLES,  0, (N-1)*(N-1)*18);
		
		// check and call events and swap buffers
		glfwSwapBuffers(window);
		glfwPollEvents(); 

		velocityStep (N, u, v, u_prev, v_prev, 100.0f, dt, vel_pos );
		densityStep (N, dens, dens_prev, u_prev, v_prev, diff, dt , vel_pos, dens_grid);
	}

	// de-allocate all resources once they've outlived their purpose
    glDeleteVertexArrays(2, vaos);
    glDeleteBuffers(2, vbos);

	// glfw: terminate, clearing all previously allocated GLFW resources.
	glfwTerminate();
	return 0;
}