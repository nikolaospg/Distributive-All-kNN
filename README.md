# Distributive-All-kNN
This is an implementation of a classical all kNN search, designed to work on a distributed memory system using OpenMPI.
In V2 the search is done using a vantage point tree, where the leaf elements are organised in buckets, and is achieved through serialisation.
The project is designed for a euclidean space and therefore the distance metric is the euclidiean distance. Maybe later a non euclidean approach could be implemented for the sake of generality.
Manipulating the readfile allows the user to be able to process his/her data in the way they like.
