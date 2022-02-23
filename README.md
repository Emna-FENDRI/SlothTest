# SlothTest

```Sloth``` is a slow-timed hash function where the resulting hash can be veriﬁed quickly.
```Sloth``` takes as input a seed and a large prime number, and is used for uncontestable random number generation, where the resulting random number is referred to as a ```beacon```.

This project consists in testing the existence of a correlation between the random value produced with a given seed and the random value produced with this same seed with one bit ﬂipped, with all other parameters unchanged.
For this, we realise a Chi-square test of independence based on sample data.
