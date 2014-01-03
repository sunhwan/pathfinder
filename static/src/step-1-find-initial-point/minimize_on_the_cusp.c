/* Starting from a point on the cusp minimize using the energy function of either 
   of the surfaces while staying on the cusp.
*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>

#include "calc_force_energy_2d.h"
#include "solve_quadratic.h"

/* Given a point, it takes the y value and solve for x on the cusp line. There will be two values
   of x. */
void Solve_for_x(double x01[2], double k_x1, double k_y1, double k_xy1, double x02[2], double k_x2, double k_y2, double k_xy2, double input_y, double *output_x_1, double *output_x_2)
{
  double a, b, c, x_1, y_1, x_2, y_2;
  int flag;

  x_1 = x01[0];
  y_1 = x01[1];
  x_2 = x02[0];
  y_2 = x02[1];

  a = k_x1 - k_x2;
  b = (-2.0*k_x1*x_1) + (k_xy1*(input_y-y_1)) + (2.0*k_x2*x_2) + (-k_xy2*(input_y-y_2));
  c = (k_x1*pow(x_1,2)) + (-k_xy1*x_1*(input_y-y_1)) + (-k_x2*pow(x_2,2)) + (k_xy2*x_2*(input_y-y_2)) + (k_y1*pow((input_y-y_1),2)) + (-k_y2*pow((input_y-y_2),2));

  flag = Solve_quadratic(a, b, c, output_x_1, output_x_2);
  if(flag != 1)
    {
      printf("Imaginary roots.\n");
      exit(1);
    }
}


int main()
{
  double x01[2], x02[2];
  double k_x1, k_y1, k_xy1, k_x2, k_y2, k_xy2;
  double tol;
  int Num_images;
  double x_cusp_ini[2], x_cusp[2], x_1[2], x_2[2];
  FILE *in, *out_e, *out, *out_p;
  int Num_iterations, iteration_number;
  double step_size;
  double potential_1, potential_2;
  double force_1[2], force_2[2];


  /* Left well surface 1.*/
  k_x1 = 2000.0;
  k_y1 = 1000.0;
  k_xy1 = -900.0;
  x01[0] = -3.0;
  x01[1] = -3.0;

  /* Right well surface 2*/
  k_x2 = 600.0;
  k_y2 = 500.0;
  k_xy2 = 500.0;
  x02[0] = 3.0;
  x02[1] = -5.0;


  in = fopen("strucutre_at_cusp", "r");
  if(in == NULL)
    {
      printf("File 'strucutre_at_cusp' not found\n");
      exit(1);
    }  
  fscanf(in, "%lf %lf ", &x_cusp_ini[0], &x_cusp_ini[1]);
  fclose(in);
  in = fopen("strucutre_at_cusp", "r");
  fscanf(in, "%lf %lf ", &x_cusp[0], &x_cusp[1]);
  fclose(in);

  /* Loop for main minimization*/
  tol = 1.0E-8;
  Num_images = 100;
  Num_iterations = 1000;
  step_size = 0.00001;
  out_e = fopen("energies_on_cusp", "w");
  out_p = fopen("points_on_cusp", "w");
  Cacl_pot_force(x_cusp, x01, k_x1, k_y1, k_xy1, &potential_1, force_1);
  Cacl_pot_force(x_cusp, x02, k_x2, k_y2, k_xy2, &potential_2, force_2);
  fprintf(out_e, "%d     %.15e  %.15e\n", 0, potential_1, potential_2);
  fprintf(out_p, "%d   %lf %lf\n", 0, x_cusp[0], x_cusp[1]);
  iteration_number = 1;
  do
    {
      /* Minimize on each surface for one step*/
      Cacl_pot_force(x_cusp, x01, k_x1, k_y1, k_xy1, &potential_1, force_1);
      x_1[0] = x_cusp[0] + (step_size * force_1[0]);
      x_1[1] = x_cusp[1] + (step_size * force_1[1]);

      Cacl_pot_force(x_cusp, x02, k_x2, k_y2, k_xy2, &potential_2, force_2); 
      x_2[0] = x_cusp[0] + (step_size * force_2[0]);
      x_2[1] = x_cusp[1] + (step_size * force_2[1]);

      printf("%lf %lf   %lf %lf\n", x_1[0], x_1[1], x_2[0], x_2[1]);     

      /* Create a new point on the cusp that is in between two points gotten from the 
         last step*/
      Find_strucutre_on_cusp_1(x01, k_x1, k_y1, k_xy1, x02, k_x2, k_y2, k_xy2, Num_images, tol, x_1, x_2, x_cusp);

      /* Energy at the new point*/
      Cacl_pot_force(x_cusp, x01, k_x1, k_y1, k_xy1, &potential_1, force_1);
      Cacl_pot_force(x_cusp, x02, k_x2, k_y2, k_xy2, &potential_2, force_2);
      fprintf(out_e, "%d     %.15e  %.15e\n", iteration_number, potential_1, potential_2);
      fprintf(out_p, "%d   %lf %lf\n", iteration_number, x_cusp[0], x_cusp[1]); 
      fflush(out_e);

      printf("%d\n", iteration_number);     

      iteration_number = iteration_number + 1;
    }while(iteration_number <= Num_iterations);
  fclose(out_e);
  fclose(out_p);

  out = fopen("minimum_on_cusp_energy_1", "w");
  fprintf(out, "%lf   %lf    %lf\n", x_cusp[0], x_cusp[1], potential_1);
  fclose(out);
  out = fopen("minimum_on_cusp_energy_2", "w");
  fprintf(out, "%lf   %lf    %lf\n", x_cusp[0], x_cusp[1], potential_2);
  fclose(out);

  return 0;
}



