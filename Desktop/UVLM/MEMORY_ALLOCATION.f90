SUBROUTINE MEMORY_ALLOCATION
USE global ! module for global variables

USE iga_global

ALLOCATE(t(N_steps))

ALLOCATE(UX(0:NR+2))

ALLOCATE(UY(0:NC+2))

ALLOCATE(GEOM(1:3,0:NDY,0:NDX))

ALLOCATE(GEOMVAL(18))

ALLOCATE(AreaElt(NR,NC))

ALLOCATE(Varea(NR*NC))

ALLOCATE(X0(NR+1,NC+1))

ALLOCATE(Y0(NR+1,NC+1))

ALLOCATE(Z0(NR+1,NC+1))

ALLOCATE(Ones(N_steps))

ALLOCATE(twist1(N_steps))

ALLOCATE(bend1(N_steps))

ALLOCATE(twist2(N_steps))

ALLOCATE(bend2(N_steps))

ALLOCATE(fore(N_steps))

ALLOCATE(Path_X(N_steps))

ALLOCATE(Path_Y(N_steps))

ALLOCATE(Path_Z(N_steps))

ALLOCATE(Roll(N_steps))

ALLOCATE(Pitch(N_steps))

ALLOCATE(XX(NR+1,NC+1,N_steps))

ALLOCATE(YY(NR+1,NC+1,N_steps))

ALLOCATE(ZZ(NR+1,NC+1,N_steps))

ALLOCATE(XR(NR+1,NC+1,N_steps))

ALLOCATE(YR(NR+1,NC+1,N_steps))

ALLOCATE(ZR(NR+1,NC+1,N_steps))

ALLOCATE(XRb(NR+1,NC+1,N_steps))

ALLOCATE(YRb(NR+1,NC+1,N_steps))

ALLOCATE(ZRb(NR+1,NC+1,N_steps))

ALLOCATE(XC(NR,NC,N_steps))

ALLOCATE(YC(NR,NC,N_steps))

ALLOCATE(ZC(NR,NC,N_steps))

ALLOCATE(XCb(NR,NC,N_steps))

ALLOCATE(YCb(NR,NC,N_steps))

ALLOCATE(ZCb(NR,NC,N_steps))

ALLOCATE(U_KIN(NR,NC,N_steps))

ALLOCATE(V_KIN(NR,NC,N_steps))

ALLOCATE(W_KIN(NR,NC,N_steps))

ALLOCATE(Outward_X(NR,NC,N_steps))

ALLOCATE(Outward_Y(NR,NC,N_steps))

ALLOCATE(Outward_Z(NR,NC,N_steps))

ALLOCATE(Tau_X(NR,NC,N_steps))

ALLOCATE(Tau_Y(NR,NC,N_steps))

ALLOCATE(Tau_Z(NR,NC,N_steps))

ALLOCATE(N_drag_X(NR,NC,N_steps))

ALLOCATE(N_drag_Y(NR,NC,N_steps))

ALLOCATE(N_drag_Z(NR,NC,N_steps))

ALLOCATE(N_lift_X(NR,NC,N_steps))

ALLOCATE(N_lift_Y(NR,NC,N_steps))

ALLOCATE(N_lift_Z(NR,NC,N_steps)) 

ALLOCATE(ALPHA(NR,NC,N_steps)) 

ALLOCATE(X_Wake(N_steps,NC+1)) 

ALLOCATE(Y_Wake(N_steps,NC+1)) 

ALLOCATE(Z_Wake(N_steps,NC+1)) 

ALLOCATE(Circ(NR,NC,N_steps))

ALLOCATE(Circ_Wake(N_steps-1,NC)) 

ALLOCATE(C_mat(NR*NC,NR*NC,N_steps))

ALLOCATE(B_mat(NR*NC,NR*NC,N_steps))

ALLOCATE(NORM_FLOW(NR,NC,N_steps))

ALLOCATE(LL(NR*NC,N_steps))

ALLOCATE(C_mat_Wake(NR*NC,(N_steps-1)*NC))

ALLOCATE(B_mat_Wake(NR*NC,(N_steps-1)*NC))

ALLOCATE(Mat(NR*NC,NR*NC))

ALLOCATE(IPIV(NEL-1))

ALLOCATE(R(NR*NC))

ALLOCATE(G(NR*NC))

ALLOCATE(VV(NEL))

ALLOCATE(MVV(NC,NR))

ALLOCATE(W_ind(NR,NC,N_steps))

ALLOCATE(VV2(NC*(N_steps-1)))

ALLOCATE(MVV2(NC,N_steps-1))

ALLOCATE(VV3(NR*NC))

ALLOCATE(MVV3(NC,NR))

ALLOCATE(Circ_foo(NR,NC,N_steps))

ALLOCATE(Circ_local(NR,NC,N_steps))

ALLOCATE(Sigma1(NR,NC,N_steps))

ALLOCATE(DFDT(NR,NC,N_steps))

ALLOCATE(DCP(NR,NC,N_steps))

ALLOCATE(POWER(N_steps))

ALLOCATE(CD(N_steps))

ALLOCATE(CL(N_steps))

ALLOCATE(C1_rollup((NC+1)*(N_steps-1),NR*NC,3))

ALLOCATE(C2_rollup((NC+1)*(N_steps-1),(N_steps-1)*NC,3))

ALLOCATE(DEL_X((NC+1)*(N_steps-1)))

ALLOCATE(DEL_Y((NC+1)*(N_steps-1)))

ALLOCATE(DEL_Z((NC+1)*(N_steps-1)))

ALLOCATE(MDEL_X((NC+1),(N_steps-1)))

ALLOCATE(MDEL_Y((NC+1),(N_steps-1)))

ALLOCATE(MDEL_Z((NC+1),(N_steps-1)))

ALLOCATE(MVV4(NC,NR))

ALLOCATE(MVV5(NC,N_steps-1))

ALLOCATE(VV4(NC*NR))

ALLOCATE(VV5(NC*(N_steps-1)))

ALLOCATE(t_spline(3*N_k))

ALLOCATE(TW_spline(3*N_k))

ALLOCATE(TW2_spline(3*N_k))

ALLOCATE(B_spline(3*N_k))

ALLOCATE(B2_spline(3*N_k))

ALLOCATE(y2it1(3*N_k))

ALLOCATE(y2it2(3*N_k))

ALLOCATE(y2ib1(3*N_k))

ALLOCATE(y2ib2(3*N_k))

ALLOCATE(Y_spline(N_steps))

ALLOCATE(Y2_spline(N_steps))

ALLOCATE(Y3_spline(N_steps))

ALLOCATE(Y4_spline(N_steps))

END SUBROUTINE MEMORY_ALLOCATION
