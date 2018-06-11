function b = createBody (Nel,Nnod);

   % -------------------------- Mesh -----------------------------
   b.Nel  = Nel;                    % Number of elements.
   b.Nnod = Nnod;                   % Number of nodes.
   b.Ndof = 6*Nnod;                 % Number of DOFs.
   b.Pn_B = zeros(b.Ndof,1);        % Undeformed nodal coordinates in the body CS.

   % ------------------------ Elements ---------------------------
   b.Lel  = zeros(Nel,1);           % Length of each element.
   b.xis  = zeros(Nel,1);           % Structural twist angle.
   b.EEs  = zeros(6,6*Nel);         % 6-by-6 matrices of cross-section stiffness.
   b.rhos = zeros(6,6*Nel);         % 6-by-6 matrices of cross-section inertia.
   b.ke_s = zeros(12,12*Nel);       % Element stiffness matrices.
   b.me_s = zeros(12,12*Nel);       % Element mass matrices.
   b.conn = zeros(3,Nel);           % Connectivity: ref node, node 1, node 2.

   % -------------------------- Other ----------------------------
   b.zeta = 0;                      % Baseline modal structural damping.
   b.Nmod = 0;                      % Number of modes to retain.
