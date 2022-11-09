density=(5) #(1 2 5 10)
thickness=(5) #(2 5 10) 


for t in "${thickness[@]}"
do
    for d in "${density[@]}"
    do

	./GEN_FOAM.py ${t} ${d}
    done
done


