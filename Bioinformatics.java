import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;
import java.util.Scanner;

import javax.print.attribute.standard.RequestingUserName;


public class Bioinformatics {
	
	public ArrayList<String> gibbsSampler(ArrayList<String> dna, int k, int t, int n){
		ArrayList<String> bestMotifs = new ArrayList<String>();
		ArrayList<String> motifs = new ArrayList<String>();
		int l = dna.get(0).length();
		for(int i = 0; i < dna.size(); i++){
			Random rn = new Random();
			int rand = rn.nextInt(l - k);
			motifs.add(dna.get(i).substring(rand, rand + k));
		}
		bestMotifs = motifs;
		for(int j = 1; j < n; j++){
			Random rn = new Random();
			int rand = rn.nextInt(t);
			ArrayList<String> temp = dna;
			temp.remove(rand);
			motifs.remove(rand);
			double[][] profile = returnProfilePseudoCounts(motifs, k);
			String motifi = ProfileMostProbableKmer(dna.get(rand), k, profile);
			motifs.add(rand, motifi);
			int score = returnScore(motifs);
			int bestscore = returnScore(bestMotifs);
			if(score < bestscore)
				bestMotifs = motifs;
		}
		return bestMotifs;
	}
	
	public ArrayList<String> randomizedMotifSearch(int k, int t, ArrayList<String> dna){
		ArrayList<String> bestMotifs = new ArrayList<String>();
		ArrayList<String> motifs = new ArrayList<String>();
		int l = dna.get(0).length();
		for(int i = 0; i < dna.size(); i++){
			Random rn = new Random();
			int rand = rn.nextInt(l - k);
			motifs.add(dna.get(i).substring(rand, rand + k));
		}
		bestMotifs = motifs;
		while(true){
			double profile[][] = returnProfilePseudoCounts(motifs, k);
			motifs = motifs(k, profile, dna);
			int score = returnScore(motifs);
			int bestscore = returnScore(bestMotifs);
			if(score < bestscore)
				bestMotifs = motifs;
			else
				return bestMotifs;
		}
	}
	
	public ArrayList<String> motifs(int k, double[][] profile, ArrayList<String> dna){
		ArrayList<String> mostProbMotifs = new ArrayList<String>();
		for(int i = 0; i < dna.size(); i++){
			String temp = ProfileMostProbableKmer(dna.get(i), k, profile);
			mostProbMotifs.add(temp);
		}
		return mostProbMotifs;
	}
	
	public ArrayList<String> greedyMotifSearch(int k, int t, ArrayList<String> dna){
		ArrayList<String> bestMotifs = new ArrayList<String>();
		for(int i = 0; i < dna.size(); i++){
			bestMotifs.add(dna.get(i).substring(0, k));
		}
		int l = dna.get(0).length();
		for(int i = 0; i <= l - k; i++){
			ArrayList<String> motifs = new ArrayList<String>();
			String temp = dna.get(0).substring(i, i + k);
			motifs.add(temp);
			for(int j = 1; j < t; j++){
				double profile[][] = returnProfile(motifs, k);
				motifs.add(ProfileMostProbableKmer(dna.get(j), k, profile));
			}
			int score = returnScore(motifs);
			int bestscore = returnScore(bestMotifs);
			if(score < bestscore)
				bestMotifs = motifs;				
		}
		return bestMotifs;
	}
	
	public ArrayList<String> greedyMotifSearchWithPseudoCounts(int k, int t, ArrayList<String> dna){
		ArrayList<String> bestMotifs = new ArrayList<String>();
		for(int i = 0; i < dna.size(); i++){
			bestMotifs.add(dna.get(i).substring(0, k));
		}
		int l = dna.get(0).length();
		for(int i = 0; i <= l - k; i++){
			ArrayList<String> motifs = new ArrayList<String>();
			String temp = dna.get(0).substring(i, i + k);
			motifs.add(temp);
			for(int j = 1; j < t; j++){
				double profile[][] = returnProfilePseudoCounts(motifs, k);
				motifs.add(ProfileMostProbableKmer(dna.get(j), k, profile));
			}
			int score = returnScore(motifs);
			int bestscore = returnScore(bestMotifs);
			if(score < bestscore)
				bestMotifs = motifs;				
		}
		return bestMotifs;
	}
	
	public String consensus(ArrayList<String> motifs){
		String cons = new String();
		char con[] = new char[motifs.get(0).length()];
		for(int i = 0; i < motifs.get(0).length(); i++){
			int count[] = new int[4];
			for(int j = 0; j < motifs.size(); j++){
				char c = motifs.get(j).charAt(i);
				switch(c){
				case 'A':
					count[0]++;
					break;
				case 'C':
					count[1]++;
					break;
				case 'G':
					count[2]++;
					break;
				case 'T':
					count[3]++;
					break;
				}
			}
			int k = 0;
			int maxIndex = 0;
			int maxCount = Integer.MIN_VALUE;
			while(k < 4){
				if(count[k] > maxCount){
					maxCount = count[k];
					maxIndex = k;
				}
				k++;
			}
			con[i] = numberToSymbol(maxIndex);
		}
		cons = String.valueOf(con);
		return cons;
	}
	
	public int returnScore(ArrayList<String> motifs){
		int score = 0;
		String cons = consensus(motifs);
		for(int i = 0; i < motifs.size(); i++){
			score += hammingDistance(cons, motifs.get(i));
		}
		return score;
	}
	
	public double[][] returnProfile(ArrayList<String> motifs, int k){
		double[][] profile = new double[4][k];
		int p = 0;
		for(int i = 0; i < motifs.size(); i++){
			String temp = motifs.get(i);
			for(int j = 0; j < temp.length(); j++){
				char c = temp.charAt(j);
				switch(c){
				case 'A':
					p = 0;
					break;
				case 'C':
					p = 1;
					break;
				case 'G':
					p = 2;
					break;
				case 'T':
					p = 3;
					break;
				}
				profile[p][j]++;
			}
		}
		for(int i = 0; i < k; i++){
			int count = 0;
			for(int j = 0; j < 4; j++){
				count += profile[j][i];
			}
			for(int j = 0; j < 4; j++){
				profile[j][i] /= count;
			}
		}
		return profile;
	}
	
	public double[][] returnProfilePseudoCounts(ArrayList<String> motifs, int k){
		double[][] profile = new double[4][k];
		int p = 0;
		for(int i = 0; i < motifs.size(); i++){
			String temp = motifs.get(i);
			for(int j = 0; j < temp.length(); j++){
				char c = temp.charAt(j);
				switch(c){
				case 'A':
					p = 0;
					break;
				case 'C':
					p = 1;
					break;
				case 'G':
					p = 2;
					break;
				case 'T':
					p = 3;
					break;
				}
				profile[p][j]++;
			}
		}
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < k; j++){
				profile[i][j]++;
			}
		}
		for(int i = 0; i < k; i++){
			int count = 0;
			for(int j = 0; j < 4; j++){
				count += profile[j][i];
			}
			for(int j = 0; j < 4; j++){
				profile[j][i] /= count;
			}
		}
		return profile;
	}
	
	public String ProfileMostProbableKmer(String text, int k, double[][] profile){
		String mostProb = new String();
		int l = text.length();
		int maxi = 0;
		double maxProb = Double.MIN_VALUE;
		double[] prob = new double[l - k + 1];
		for(int i = 0; i <= l - k; i++){
			String temp = text.substring(i, i + k);
			prob[i] = calcProb(temp, profile);
		}
		for(int i = 0; i < prob.length; i++){
			if(prob[i] > maxProb){
				maxProb = prob[i];
				maxi = i;
			}
		}
		mostProb = text.substring(maxi, maxi + k);
		return mostProb;
	}
	
	public double calcProb(String temp, double[][] profile){
		double ans = 1.0;
		int p = 0;
		for(int i = 0; i < temp.length(); i++){
			char c = temp.charAt(i);
			switch(c){
			case 'A':
				p = 0;
				break;
			case 'C':
				p = 1;
				break;
			case 'G':
				p = 2;
				break;
			case 'T':
				p = 3;
				break;
			}
			ans *= profile[p][i];
		}
		return ans;
	}
	
	public String medianString(ArrayList<String> dna, int k){
		String median = new String();
		int distance = Integer.MAX_VALUE;
		for(int i = 0; i < Math.pow(4, k); i++){
			String pattern = numberToPattern(i, k);
			if(distance > distanceBetweenPatternAndStrings(pattern, dna)){
				distance = distanceBetweenPatternAndStrings(pattern, dna);
				median = pattern;
			}
		}
		return median;
	}
	
	public int distanceBetweenPatternAndStrings(String pattern, ArrayList<String> dna){
		int distance = 0;
		int k = pattern.length();
		for(int i = 0; i < dna.size(); i++){
			int hamDist = Integer.MAX_VALUE;
			for(int j = 0; j <= dna.get(i).length() - k; j++){
				if(hamDist > hammingDistance(pattern, dna.get(i).substring(j, j + k)))
					hamDist = hammingDistance(pattern, dna.get(i).substring(j, j + k));
			}
			distance += hamDist;
		}
		return distance;
	}
	
	public ArrayList<String> motifEnumeration(ArrayList<String> dna, int k, int d){
		ArrayList<String> motif = new ArrayList<String>();
		ArrayList<String> kmer = new ArrayList<String>();
		int l = dna.size();
		for(int i = 0; i < l; i++){
			String temp = dna.get(i);
			for(int n = 0; n <= temp.length() - k; n++)
				kmer.add(temp.substring(n, n + k));
		}
		for(int i = 0; i < kmer.size(); i++){
			ArrayList<String> p = neighbors(kmer.get(i), d);
			for(int j = 0; j < p.size(); j++){
				int count = 0;
				for(int m = 0; m < l; m++){
					for(int n = 0; n <= dna.get(m).length() - k; n++){
						ArrayList<String> t = neighbors(dna.get(m).substring(n, n + k),d);
						if(t.contains(p.get(j))){
							count++;
							break;
						}
					}
				}
				if(count == l)
					motif.add(p.get(j));
			}
		}
		for(int i = 0; i < motif.size(); i++){
			for(int j = 0; j < motif.size(); j++){
				if(motif.get(i).equalsIgnoreCase(motif.get(j)) && i != j){
					motif.remove(j);
					j--;
				}
			}
		}
		return motif;
	}
	
	public ArrayList<String> frequentWordsWithMismatches(String text, int k ,int d){
		ArrayList<String> freq = new ArrayList<String>();
		int l = text.length();
		int maxCount = 0;
		int[] close = new int[(int) Math.pow(4, k)];
		int[] freqarray = new int[(int) Math.pow(4, k)];
		for(int i = 0; i <= l - k; i++){
			ArrayList<String> neig = neighbors(text.substring(i, i + k), d);
			for(int j = 0; j < neig.size(); j++){
				int index = patternToNumber(neig.get(j));
				close[index] = 1;
			}
		}
		for(int i = 0; i < Math.pow(4, k); i++){
			if(close[i] == 1){
				String pattern = numberToPattern(i, k);
				freqarray[i] = approxPatternCount(text, pattern, d);
			}	
		}
		for(int i = 0; i < freqarray.length; i++){
			if(maxCount < freqarray[i])
				maxCount = freqarray[i];
		}
		for(int i = 0; i < Math.pow(4, k); i++){
			if(freqarray[i] == maxCount){
				String pattern = numberToPattern(i, k);
				freq.add(pattern);
			}	
		}
		return freq;
	}
	
	public ArrayList<String> neighbors(String pattern, int d){
		ArrayList<String> neig = new ArrayList<String>();
		int p = pattern.length();
		ArrayList<String> nucl = new ArrayList<String>(Arrays.asList("A", "C", "G", "T"));
		if(d == 0){
			neig.add(pattern);
			return neig;
		}
		if(pattern.length() == 1){
			return nucl;
		}
		String suffPatt = pattern.substring(1, p);
		ArrayList<String> suffixNeig = neighbors(suffPatt, d);
		for(int i = 0; i < suffixNeig.size(); i++){
			String temp = suffixNeig.get(i);
			if(hammingDistance(suffPatt, temp) < d){
				for(int j = 0; j < 4; j++)
					neig.add(nucl.get(j)+ temp);
			}
			else
				neig.add((pattern.substring(0, 1)) + temp);
		}
		return neig;
	}
	
	public ArrayList<Integer> minSkew(String text){
		ArrayList<Integer> min = new ArrayList<Integer>();
		int[] sk = skew(text);
		int mins = Integer.MAX_VALUE;
		for(int i = 0; i < sk.length; i++){
			if(mins > sk[i]){
				mins = sk[i];
			}
		}
		for(int i = 0; i < sk.length; i++){
			if(sk[i] == mins)
				min.add(i);
		}
		return min;
	}
	
	public int hammingDistance(String a, String b){
		int l = a.length();
		int ham = 0;
		for(int i = 0; i < l; i++){
			if(a.charAt(i) != b.charAt(i))
				ham++;
		}
		return ham;
	}
	
	public int[] computeFreq(String text, int k){
		int[] freq = new int[(int) Math.pow(4, k)];
		for(int i = 0; i <= text.length() - k; i++){
			String temp = text.substring(i, i + k);
			int j = patternToNumber(temp);
			freq[j] += 1; 
		}
		
		return freq;
	}
	
	public int[] skew(String text){
		int l = text.length();
		int[] sk = new int[l+1];
		sk[0] = 0;
		for(int i = 0; i < l; i++){
			char c = text.charAt(i);
			if(c == 'C')
				sk[i+1] = sk[i] - 1;
			else if (c == 'G')
				sk[i+1] = sk[i] + 1;
			else
				sk[i+1] = sk[i];
		}
		return sk;
	}
	
	public int patternToNumber(String text){
		int l = text.length();
		if(l == 0)
			return 0;
		String prefix = text.substring(0, l - 1);
		char symbol = text.charAt(l - 1);
		return 4 * patternToNumber(prefix) + symbolToNumber(symbol);		
	}
	
	public String numberToPattern(int index, int k){
		if(k == 1)
			return String.valueOf(numberToSymbol(index));
		int prefixIndex = index / 4;
		int r = index % 4;
		char symbol = numberToSymbol(r);
		String prefixPattern = numberToPattern(prefixIndex, k - 1);
		return prefixPattern + symbol;
	}
	
	public int symbolToNumber(char c){
		switch(c){
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		}
		return -1;
	}
	
	public char numberToSymbol(int c){
		switch(c){
		case 0:
			return 'A';
		case 1:
			return 'C';
		case 2:
			return 'G';
		case 3:
			return 'T';
		}
		return 'Z';
	}
		
	public ArrayList<Integer> patternMatching(String text, String pattern){
		int p = pattern.length();
		int t = text.length();
		ArrayList<Integer> ans = new ArrayList<Integer>();
		for(int i = 0; i <= t - p; i++){
			String temp = text.substring(i, i + p);
			if(temp.equalsIgnoreCase(pattern))
				ans.add(i);
		}
		return ans;
	}
	
	public ArrayList<Integer> approxPatternMatching(String  pattern, String text, int d){
		ArrayList<Integer> ans =  new ArrayList<Integer>();
		int l = text.length();
		int p = pattern.length();
		for(int i = 0; i <= l - p; i++){
			String temp = text.substring(i, i + p);
			if(hammingDistance(pattern, temp) <= d)
				ans.add(i);
		}
		return ans;
	}
	
	public int approxPatternCount(String text, String pattern, int d){
		int count = 0;
		int l = text.length();
		int p = pattern.length();
		for(int i = 0; i <= l - p; i++){
			String temp = text.substring(i, i + p);
			if(hammingDistance(pattern, temp) <= d)
				count++;
		}
		return count;
	}
	
	public int patternCount(String s, String p) {
		int count = 0;
		int c = s.length()- p.length();
		int modp = p.length();
		for(int i = 0; i <= c; i++){
			String t = s.substring(i, i + modp);
			if(t.equalsIgnoreCase(p))
				count++;
		}
		return count;
	}
	
	public String reverseComplement(String text){
		int l = text.length();
		char[] rc = new char[l];
		int j = 0;
		for(int i = l-1; i >= 0; i--){
			char c = text.charAt(i);
			switch(c){
				case 'A':
					rc[j] = 'T';
					break;
				case 'T':
					rc[j] = 'A';
					break;
				case 'C':
					rc[j] = 'G';
					break;
				case 'G':
					rc[j] = 'C';
					break;
			}
			j++;
		}
		String ans = new String(rc);
		return ans;
	}
	
	public ArrayList<String> frequentPatterns(String text, int k){
		ArrayList<String> fp = new ArrayList<String>();
		int c = text.length();
		int[] count = new int[c];
		int maxCount = 0;
		for(int i = 0; i <= c - k; i++) {
			String temp = text.substring(i, i + k);
			count[i] = patternCount(text, temp);
			if(count[i] > maxCount)
				maxCount = count[i];
		}
		for(int i = 0; i <= c - k; i++) {
			if(count[i] == maxCount)
				fp.add(text.substring(i, i + k));
		}
		for(int i = 0; i < fp.size(); i++){
			for(int j = 0; j < fp.size(); j++){
				if(fp.get(i).equalsIgnoreCase(fp.get(j)) && i != j)
					fp.remove(j);
			}
		}
		return fp;
	}
	
	public ArrayList<String> findClump(String genome, int k, int l, int t){
		ArrayList<String> fc = new ArrayList<String>();
		int c = genome.length();
		for(int i = 0; i <= c - l; i++) {
			int[] count = new int[l-k+1];
			String window = genome.substring(i, i + l);
			for(int j = 0; j <= l - k; j++){
				String temp = genome.substring(i + j, i + j + k);
				count[j] = patternCount(window, temp);
				if(count[j] >= t && !fc.contains(temp))
					fc.add(temp);
			}
		}
		return fc;
	}

	public ArrayList<String> frequentWordsWithMismatchesAndReverse(String text, int k ,int d){
		ArrayList<String> freq = new ArrayList<String>();
		int l = text.length();
		int maxCount = 0;
		int[] close = new int[(int) Math.pow(4, k)];
		int[] freqarray = new int[(int) Math.pow(4, k)];
		for(int i = 0; i <= l - k; i++){
			ArrayList<String> neig = neighbors(text.substring(i, i + k), d);
			for(int j = 0; j < neig.size(); j++){
				int index = patternToNumber(neig.get(j));
				close[index] = 1;
			}
		}
		for(int i = 0; i < Math.pow(4, k); i++){
			if(close[i] == 1){
				String pattern = numberToPattern(i, k);
				String rev = reverseComplement(pattern);
				freqarray[i] = approxPatternCount(text, pattern, d) + approxPatternCount(text, rev, d);
			}	
		}
		for(int i = 0; i < freqarray.length; i++){
			if(maxCount < freqarray[i])
				maxCount = freqarray[i];
		}
		for(int i = 0; i < Math.pow(4, k); i++){
			if(freqarray[i] == maxCount){
				String pattern = numberToPattern(i, k);
				freq.add(pattern);
			}	
		}
		return freq;
	}
}

