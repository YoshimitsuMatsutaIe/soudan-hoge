Write-Output "コンパイル実行中..."

Set-Location .\o\soft_robot\derivation_of_kinematics

Write-Output "Phiのコンパイル実行中..."
gcc -shared -fPIC -o ./derived/so/Phi_0.so ./derived/c_src/Phi_s/Phi_0.c
gcc -shared -fPIC -o ./derived/so/Phi_1.so ./derived/c_src/Phi_s/Phi_1.c
gcc -shared -fPIC -o ./derived/so/Phi_2.so ./derived/c_src/Phi_s/Phi_2.c
gcc -shared -fPIC -o ./derived/so/Phi_3.so ./derived/c_src/Phi_s/Phi_3.c
gcc -shared -fPIC -o ./derived/so/Phi_4.so ./derived/c_src/Phi_s/Phi_4.c


Write-Output "Jのコンパイル実行中..."
gcc -shared -fPIC -o ./derived/so/J_0.so ./derived/c_src/J_s/J_0.c
gcc -shared -fPIC -o ./derived/so/J_1.so ./derived/c_src/J_s/J_1.c
gcc -shared -fPIC -o ./derived/so/J_2.so ./derived/c_src/J_s/J_2.c
gcc -shared -fPIC -o ./derived/so/J_3.so ./derived/c_src/J_s/J_3.c
gcc -shared -fPIC -o ./derived/so/J_4.so ./derived/c_src/J_s/J_4.c


Write-Output "完了!"