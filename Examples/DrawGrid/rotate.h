/*********************************
 Implementation of quaternion-based
 object rotation
 
 Functions
 clickmotion - call when holded mouse moved
 click       - call when user clicks
 motion      - call when mouse just moves
 quatinit    - flush rotation value
 rotate      - multiply GL matrix

 Dependency: rotate.h,
   Standard: math.h 
   Specific: glut.h
**********************************/

#ifndef _ROTATE_H
#define _ROTATE_H

void clickmotion(int nmx, int nmy);
void motion(int nmx, int nmy);
void click(int b, int s, int nmx, int nmy);
void quatget(double * vec);
void quatinit();
void quatpush();
void quatpop();
void rotate();
void rotate_from_stack();
void rotatevector(double * vec);
void reverse_rotatevector(double * vec);
void revrotatevector(double * vec);
void rotatevector_from_stack(double * vec);
void reverse_rotatevector_from_stack(double * vec);

void quatget4(double gq[4]);
void quatget4_from_stack(double gq[4]);
void quatset4(double gq[4]);

#endif
