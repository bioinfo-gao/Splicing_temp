cat ~/.ssh/id_rsa.pub | ssh user@123.45.56.78 "mkdir -p ~/.ssh && cat >>  ~/.ssh/authorized_keys"
cat ~/.ssh/id_rsa.pub | ssh zgao1@edge.hpc.biogen.com "mkdir -p ~/.ssh && cat >>  ~/.ssh/authorized_keys"
cat ~/.ssh/id_rsa.pub | ssh zgao1@camhpcve23.hpc.biogen.com "mkdir -p ~/.ssh && cat >>  ~/.ssh/authorized_keys"