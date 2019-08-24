struct rb_state {
  vector v;
  vector w;
  vector P;
  vector4 O; // orientation
}

rb_state rbd_sate_update
(const float dt;
 const vector dv;
 const vector dw;
 const vector dP;
 const vector dO;
 const rb_state rb0)
{
  rb_state rb1;
  rb1.v = rb0.v + dt*dv;
  rb1.w = rb0.w + dw*dt;
  rb1.P = rb0.P + dt*dP;
  rb1.O = qmultiply(rb0.O, quaternion(dt,dO) );
  return rb1;
}

void DampingForce
(vector F;
 vector T;
 const float rho;
 const vector V;
 const vector W)
{
  float lenW = length(W);
  T += (-0.5*rho*lenW*5.0)*W;
  float lenV = length(V);
  if( lenV < 1.e-10 ){ return; }
  vector uV = V/lenV;
  F += -0.5*rho*uV*lenV*lenV;
}

void velocity_rigid_body
(vector dv;
 vector dw;
 const float mass;
 const matrix3 I;
 const vector gravity;
 const float rho;
 const rb_state rb)
{
  vector4 Oinv = qinvert(rb.O);
  vector V = qrotate( Oinv, rb.v );
  vector F = qrotate( Oinv, gravity)*mass;
  vector T = set(0,0,0);
  DampingForce(F,T,
               rho,
               V,rb.w);
  dv = qrotate( rb.O, F ) / mass;
  dw = invert(I)*(cross(rb.w, I*rb.w)+T);
}


void rbd_vacuum_forwardeular
(vector v;
 vector w;
 vector P;
 vector4 orient;
 const float dt;
 const float mass;
 const matrix3 I;
 const vector gravity;
 const float rho)
{
  rb_state rb0 = rb_state(v,w,P,orient);
  ////
  vector dv,dw;
  velocity_rigid_body(dv,dw,
                      mass,I,gravity,rho, rb0);
  rb_state rb1 = rbd_sate_update(dt, dv,dw,rb0.v,rb0.w, rb0);
  v = rb1.v;
  w = rb1.w;
  P = rb1.P;
  orient = rb1.O;
}


void rbd_vacuum_rungekutta
(vector v;
 vector w;
 vector P;
 vector4 orient;
 const float dt;
 const float mass;
 const matrix3 I;
 const vector gravity;
 const float rho)
{
  rb_state rb0 = rb_state(v,w,P,orient);
  ////
  vector dv1,dw1;
  velocity_rigid_body(dv1,dw1,
                      mass,I,gravity,rho,   rb0);
  rb_state rb1 = rbd_sate_update(dt*0.5, dv1,dw1,rb0.v,rb0.w, rb0);
  ////
  vector dv2,dw2;
  velocity_rigid_body(dv2,dw2,
                      mass,I,gravity,rho,   rb1);
  rb_state rb2 = rbd_sate_update(dt*0.5, dv2,dw2,rb1.v,rb1.w, rb0);
  ////
  vector dv3,dw3;
  velocity_rigid_body(dv3,dw3,
                      mass,I,gravity,rho,   rb2);
  rb_state rb3 = rbd_sate_update(dt*1.0, dv3,dw3,rb2.v,rb2.w, rb0);
  ////
  vector dv4,dw4;
  velocity_rigid_body(dv4,dw4,
                      mass,I,gravity,rho,   rb3);
  ////
  vector dv = (dv1 + 2*dv2 + 2*dv3 + dv4)/6.0;
  vector dw = (dw1 + 2*dw2 + 2*dw3 + dw4)/6.0;
  vector dP = (rb0.v + 2*rb1.v + 2*rb2.v + rb3.v)/6.0;
  vector dO = (rb0.w + 2*rb1.w + 2*rb2.w + rb3.w)/6.0;
  rb_state rb5 = rbd_sate_update(dt, dv,dw,dP,dO, rb0);
  v = rb5.v;
  w = rb5.w;
  P = rb5.P;
  orient = rb5.O;
}
