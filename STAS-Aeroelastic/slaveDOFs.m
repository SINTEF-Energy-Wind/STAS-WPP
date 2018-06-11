function slv = slaveDOFs (idofs)

slv = [idofs(2)+[1:6] idofs(3)+[1:6] idofs(4)+[1:6] idofs(5)+[1:3] ...
       idofs(6)+[1:6] idofs(7)+[1:6] idofs(8)+[1:6]].';