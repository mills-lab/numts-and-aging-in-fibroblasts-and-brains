#INPUT: 1KG.txt, merged.vcf, sample.list, genome.list, & *.list
#REQUIRE: inter.0118.o

cat ../data/merged.vcf | grep "#" > header
cat ../data/merged.vcf | grep  -v "#" | grep "PASS" > PASS.mergedAll.vcf
less  PASS.mergedAll.vcf | awk '{print $1,$2,$2,$0}' > Q.txt
less ../data/1KG.txt | awk '{print $1,$2,$2,$0}' > S.txt

cp ../data/inter.0118.o .
./inter.0118.o

less inter.txt | awk '{if($12=="") print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > filtered.PASS.mergedAll
cat header filtered.PASS.mergedAll > filtered.PASS.mergedAll.vcf
less inter.txt | awk '{if($12!="") print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > 1KG.overlapped.PASS.mergedAll
cat header 1KG.overlapped.PASS.mergedAll > 1KG.overlapped.PASS.mergedAll.vcf


while read -r genome
do
cat ../data/${genome}.list ../data/sample.list | sort -n | uniq -c | awk '{if($1==1) print $2}' > non-overlap.list

rm ${genome}.list.process.1
while read -r a 
do
less filtered.PASS.mergedAll.vcf | grep ${a} >> ${genome}.list.process.1
done < ../data/${genome}.list

less ${genome}.list.process.1| sort -n | uniq > ${genome}.list.process.2

rm non-overlap.list.process.1
while read -r a
do 
less filtered.PASS.mergedAll.vcf | grep ${a} >> non-overlap.list.process.1
done < non-overlap.list

less non-overlap.list.process.1 | sort -n | uniq > non-overlap.list.process.2


cat ${genome}.list.process.2  non-overlap.list.process.2| sort -n | uniq -c | awk '{if($1==1) print $0}' > ${genome}.list.process.3

rm ${genome}.list.process.4
while read -r a 
do
less ${genome}.list.process.3 | grep ${a} >> ${genome}.list.process.4
done < ../data/${genome}.list

less ${genome}.list.process.4 | sort -n | uniq | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}'> ${genome}.list.specific


##OUTPUT1
cat header ${genome}.list.specific > ${genome}.list.specific.vcf

while read -r a 
do
less ${genome}.list.process.3 | grep ${a} | sort -n | uniq | awk '{print $2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9}' > ${genome}.list.${a}.specific


##OUTPUT2
cat header ${genome}.list.${a}.specific > ${genome}.list.${a}.specific.vcf

rm non-overlap.list

done < ../data/${genome}.list

done < ../data/genome.list


cat *.list.specific | awk '{print $1,$2,$2,$0}' > S.txt

./inter.0118.o

less inter.txt | awk '{if($12=="") print $4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11}' > non-genome-specific.PASS.mergedAll


##OUTPUT3
cat header  non-genome-specific.PASS.mergedAll >  non-genome-specific.PASS.mergedAll.vcf