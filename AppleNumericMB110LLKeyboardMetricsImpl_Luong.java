/**
 * @author Luong
 * Professor Kart
 * COSC 3327
 * Challenge 7
 * 12-01-2019
 */
package keyboard;

import static keyboard.Key.*;
import static keyboard.KeyLayout.COLEMAK;
import static keyboard.KeyLayout.DVORAK;
import static keyboard.KeyLayout.QWERTY;
import static keyboard.KeyLayout.ROTATION_13;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import combinatorics.Permutation;


public class AppleNumericMB110LLKeyboardMetricsImpl_Luong 
							implements KeyboardMetrics {
	private List<Key> vertexLabels;
	private int[][] adjacencyMatrix;
	private int[][] distanceMatrix;
	private Key homeKey;
	
	private static Map<KeyLayout, Key> keyLayoutToHomeKeyMap;
	private static Map<KeyLayout, Map<Key, Set<Key>>> keyLayoutToKeyToNeighborMapMap;
	
	static
	{
		keyLayoutToHomeKeyMap = new HashMap<KeyLayout, Key>();
		keyLayoutToHomeKeyMap.put(QWERTY, J);
		keyLayoutToHomeKeyMap.put(DVORAK, H);
		keyLayoutToHomeKeyMap.put(COLEMAK, N);
		keyLayoutToHomeKeyMap.put(ROTATION_13, W);
		
		keyLayoutToKeyToNeighborMapMap = new HashMap<KeyLayout, Map<Key, Set<Key>>>();
		Map<Key, Set<Key>> keyToNeighborMap_QWERTY = getKeyToNeighborMap_QWERTY();
//		Map<Key, Set<Key>> keyToNeighborMap_DVORAK = applyPermutationToMap(keyToNeighborMap_QWERTY, getQWERTYToDvorakPermutation());
//		Map<Key, Set<Key>> keyToNeighborMap_COLEMAK = applyPermutationToMap(keyToNeighborMap_QWERTY, getQWERTYToColemakPermutation());
//		Map<Key, Set<Key>> keyToNeighborMap_ROT_13 = applyPermutationToMap(keyToNeighborMap_QWERTY, getQWERTYToRotation13Permutation());
		keyLayoutToKeyToNeighborMapMap.put(QWERTY, keyToNeighborMap_QWERTY);
//		keyLayoutToKeyToNeighborMapMap.put(DVORAK, keyToNeighborMap_DVORAK);
//		keyLayoutToKeyToNeighborMapMap.put(COLEMAK, keyToNeighborMap_COLEMAK);
//		keyLayoutToKeyToNeighborMapMap.put(ROTATION_13, keyToNeighborMap_ROT_13);
	}
	
	//Constructor:
	public AppleNumericMB110LLKeyboardMetricsImpl_Luong(KeyLayout keyLayout)
	{
		this.homeKey = keyLayoutToHomeKeyMap.get(keyLayout);
		Map<Key, Set<Key>> keyToNeighborsMap = keyLayoutToKeyToNeighborMapMap.get(keyLayout);
		init(keyToNeighborsMap, new ArrayList<Key>(keyToNeighborsMap.keySet()));
	}
	
	public void init(Map<Key, Set<Key>> physicalKeyToNeighborsMap, List<Key> vertexLabels)
	{
		this.vertexLabels = vertexLabels;
		this.adjacencyMatrix = getAdjacencyMatrix(physicalKeyToNeighborsMap, vertexLabels);
		this.distanceMatrix = getDistanceMatrix(adjacencyMatrix);
	}
	
	//Helper method: Building the adjacency matrix for the keyboard.
	//Purpose: finding the neighbors of each keys.
	private static int[][] getAdjacencyMatrix(Map<Key, Set<Key>> physicalKeyToNeighborsMap, List<Key> vertexLabels)
	{
		assert physicalKeyToNeighborsMap.keySet().equals(new HashSet<Key>(vertexLabels)) : "vertexLabels inconsistent with physicalKeyToNeighborsMap! : vertexLabels = " + vertexLabels + " physicalKeyToNeighborsMap.keySet() = " + physicalKeyToNeighborsMap.keySet();
		final int SIZE = physicalKeyToNeighborsMap.keySet().size();
		int[][] adjacencyMatrix = new int[SIZE][SIZE];
		
		//Adjacency matrix building:
		
		for(int i = 0; i < SIZE; i++)	//1st pointer to check for the column values in the matrix
		{
			for(int j = 0; j< SIZE; j++)	//2nd pointer to check for the row values in the matrix
			{
				if(physicalKeyToNeighborsMap.get(vertexLabels.get(i)).contains(vertexLabels.get(j))) // if the two pointers lands at the same vertex.
				{
					//there is an edge-relationship between the two nodes
					adjacencyMatrix[i][j] = 1;
				}
				
				else
				{
					//No relationship found!
					adjacencyMatrix[i][j] = 0;
				}
			}
		}//end for!
		
		return adjacencyMatrix;
	}//end of getAdjacenctMatrix().
	
	//Matrix multiplication
	private static int[][] multiply(int[][] A, int[][] B)
	{
		int rowCount_A = A.length;
		assert rowCount_A > 0 : "rowCount_A = 0!";
		int columnCount_A = A[0].length;
		int rowCount_B = B.length;
		assert rowCount_B > 0 : "rowCount_B = 0!";
		int columnCount_B = B[0].length;
		assert columnCount_A == rowCount_B : "columnCount_A = " + columnCount_A + " <> " + rowCount_B + " = rowCount_B!";
		
		int[][] C = new int[rowCount_A][columnCount_B];
        for (int i = 0; i < rowCount_A; i++)
            for (int j = 0; j < columnCount_B; j++)
                for (int k = 0; k < columnCount_A; k++)
                    C[i][j] += A[i][k] * B[k][j];
        
        return C;
	}//end of multiply().
	
	//Helper method: The distance between each keys matrix
	//This is found based on the adjacencyMatrix.
	private static int[][] getDistanceMatrix(int[][] adjacencyMatrix)
	{
		int vertexCount = adjacencyMatrix.length;
		assert vertexCount > 0 : "rowCount = 0!";
		int[][] distanceMatrix = new int[vertexCount][vertexCount];
		
		//making a copy of adjacencyMatrix to manipulate its values. 
		int[][] copiedAdjacency = adjacencyMatrix;
		
		
		//Looping through the distance with the end values of 14 
		//as there is only 14 possible distance between two keys (as of now!).
		//Might not guarantee for degenerate cases. 
		for(int index = 1; index <= 14; index++)
		{
			for(int i = 0; i < vertexCount; i++)	//Looping through the columns.
			{
				for(int j = 0; j < vertexCount; j++)	//looping through the rows.
				{
					//if two pointer do not point to the same [row][column]
					//while adjacency matrix is occupied at i and j and the desired distance matrix is not.
					if(i!=j && copiedAdjacency[i][j] != 0 && distanceMatrix[i][j] == 0) 
					{
						distanceMatrix[i][j] = index;
					}
				}
			}
			
			//At each iteration of index, the adjacency matrix multiply by itself 
			//and these values will determine the distance between each keys.
			copiedAdjacency = multiply(copiedAdjacency,adjacencyMatrix);
			
		}
		return distanceMatrix;
	}//end of getDistanceMatrix().
	
	/* (non-Javadoc)
	 * @see keyboard.KeyboardMeasurements#getDistance(keyboard.PhysicalKey, keyboard.PhysicalKey)
	 */
	@Override
	//purpose: finding the distance between two specific key.
	public double getDistance(Key key1, Key key2) {
		int index1 = getIndex(vertexLabels, key1);
		int index2 = getIndex(vertexLabels, key2);
		return distanceMatrix[index1][index2];
	}//end of getDistance().

	private static <E> int getIndex(List<E> list, E element)
	{
		boolean foundIndex = false;
		int i = 0;
		while(!foundIndex && i < list.size())
		{
			foundIndex = (list.get(i) == element);
			if(!foundIndex) i++;
		}
		int rv = -1;
		if(foundIndex) rv = i;
		
		return rv;
	}

	@Override
	//purpose: this method will return the distance between each keys, when they are passed in as a string form.
	public double getDistance(String str) {
		
		double distance = 0;
		Key currentKey = homeKey;
		
		//looping through the string of characters.
		for(int i = 0; i < str.length(); i++)
		{
			//Looping through each character in a string
			Set<Key> currentKeySet = getKeySet(str.charAt(i));
			
			//Finding the closest neighbor between each key...
			Key nextDoorKey = getClosestKey(currentKeySet, currentKey);
			//as each distance between each key is found, we keep on adding them until the last 
			//character in string.
			distance += getDistance(currentKey, nextDoorKey);
			currentKey = nextDoorKey;
		}
		
		
		return distance;
	}//end of getDistance(String str).
	
	//Helper method: finding the closest key next to the key being assigned.
	//Essentially a swap function. 
	private Key getClosestKey(Set<Key> keySet, Key key)
	{
		
		List<Key> keyList = new ArrayList<Key>(keySet);
		double minDistance = getDistance(keyList.get(0), key);

		Key minDistanceKey = null;
		
		//looping through the key layouts
		for(int i = 0; i < keyList.size(); i++)
		{
			//Finding the distance between each key in the set
			double keyDistance = getDistance(keyList.get(i), key);
			//swap place with the smallest values.
			if(keyDistance <= minDistance)
			{
				//swaps:
				minDistance = keyDistance;
				
				minDistanceKey = keyList.get(i);
				
			}//end if.
		}//end for.
		
		
		return minDistanceKey;
	}//end of getClosestKey().

	private static Set<Key> getKeySet(char character)
	{
		List<Key> keyList = Arrays.asList(Key.values());
		Set<Key> characterProducingKeysSet = new HashSet<Key>();
		for(int i = 0; i < keyList.size(); i++)
		{
			Key key = keyList.get(i);
			assert key != null : "key is null!";
			boolean keyProducesCharacter = (key.getNormalCharacter() != null && key.getNormalCharacter() == character) || (key.getShiftModifiedCharacter() != null && key.getShiftModifiedCharacter() == character);
			if(keyProducesCharacter) characterProducingKeysSet.add(key);
		}
		return characterProducingKeysSet;
	}//end of getKeySet().
	
	//helper method:
	//This will create a map that put a particular key with its neighbors into 
	// a map or dictionary.
	private static Map<Key, Set<Key>> getKeyToNeighborMap_QWERTY()
	{
		Map<Key, Set<Key>> keyToNeighborSetMap = new HashMap<Key, Set<Key>>();
		
		//Produce keyToNeighborSetMap here
		//You might want to take a look at getSet()
		//CHARACTER KEYS:
		keyToNeighborSetMap.put(A, getSet(SHIFT_1, Z, S, W, Q));
		keyToNeighborSetMap.put(S, getSet(A,Z,X,D,E,W));
		keyToNeighborSetMap.put(D, getSet(S,X,C,F,R,E));
		keyToNeighborSetMap.put(F, getSet(D,C,V,G,T,R));
		keyToNeighborSetMap.put(G, getSet(F,V,B,H,Y,T));
		keyToNeighborSetMap.put(H, getSet(Y,G,B,N,J,U));
		keyToNeighborSetMap.put(J, getSet(H, N, M, K,I,U));
		keyToNeighborSetMap.put(K, getSet(I, J, M, COMMA, L, O));
		keyToNeighborSetMap.put(L, getSet(O, K, COMMA, PERIOD, SEMICOLON, P));
		keyToNeighborSetMap.put(Z, getSet(SHIFT_1, A, S, X));
		keyToNeighborSetMap.put(X, getSet(Z,S,D,C));
		keyToNeighborSetMap.put(C, getSet(X,D,F,V,SPACEBAR_1));
		keyToNeighborSetMap.put(V, getSet(C,F,G,B,SPACEBAR_2));
		keyToNeighborSetMap.put(B, getSet(V,G,H,N,SPACEBAR_3));
		keyToNeighborSetMap.put(N, getSet(B,H,J,M,SPACEBAR_4));
		keyToNeighborSetMap.put(M, getSet(N,J,K,COMMA,SPACEBAR_5));
		keyToNeighborSetMap.put(Q, getSet(TAB, ONE, TWO, W, A));
		keyToNeighborSetMap.put(W, getSet(TWO,THREE,E,S,A,Q));
		keyToNeighborSetMap.put(E, getSet(THREE,FOUR,R,D,S,W));
		keyToNeighborSetMap.put(R, getSet(FOUR,FIVE,T,F,D,E));
		keyToNeighborSetMap.put(T, getSet(FIVE, R, F, G, Y, SIX));
		keyToNeighborSetMap.put(Y, getSet(SIX, T, G, H, U, SEVEN));
		keyToNeighborSetMap.put(U, getSet(SEVEN, Y, H, J, I, EIGHT));
		keyToNeighborSetMap.put(I, getSet(EIGHT, U, J, K, O, NINE));
		keyToNeighborSetMap.put(O, getSet(NINE, I, K, L, P, ZERO));
		keyToNeighborSetMap.put(P, getSet(ZERO, O, L, SEMICOLON, LEFT_BRACKET, MINUS));
		

		//FUNCTIONAL KEYS:
		keyToNeighborSetMap.put(SHIFT_1, getSet(Z,A));
		keyToNeighborSetMap.put(SHIFT_2, getSet(RETURN, TICK, FORESLASH));
		keyToNeighborSetMap.put(BACKTICK, getSet(TAB,ONE));
		keyToNeighborSetMap.put(TAB, getSet(ONE, Q, BACKTICK));
		keyToNeighborSetMap.put(SPACEBAR_1, getSet(C));
		keyToNeighborSetMap.put(SPACEBAR_2, getSet(V));
		keyToNeighborSetMap.put(SPACEBAR_3, getSet(B));
		keyToNeighborSetMap.put(SPACEBAR_4, getSet(N));
		keyToNeighborSetMap.put(SPACEBAR_5, getSet(M));
		keyToNeighborSetMap.put(BACKSLASH, getSet(RIGHT_BRACKET, RETURN));
		keyToNeighborSetMap.put(RETURN, getSet(BACKSLASH, RIGHT_BRACKET, TICK, SHIFT_2));
		keyToNeighborSetMap.put(TICK, getSet(RETURN, SHIFT_2, FORESLASH, SEMICOLON, LEFT_BRACKET, RIGHT_BRACKET));
		keyToNeighborSetMap.put(RIGHT_BRACKET, getSet(EQUALS, LEFT_BRACKET, TICK, RETURN, BACKSLASH));
		keyToNeighborSetMap.put(LEFT_BRACKET, getSet(RIGHT_BRACKET, TICK, SEMICOLON, P, MINUS, EQUALS));
		keyToNeighborSetMap.put(SEMICOLON, getSet(LEFT_BRACKET, TICK, FORESLASH, PERIOD, L, P));
		keyToNeighborSetMap.put(FORESLASH, getSet(PERIOD, SEMICOLON, TICK, SHIFT_2));
		keyToNeighborSetMap.put(PERIOD, getSet(COMMA, L, SEMICOLON, FORESLASH));
		keyToNeighborSetMap.put(COMMA, getSet(M,K,L,PERIOD));
		
		//Numbers + operators key: 
		keyToNeighborSetMap.put(ONE, getSet(BACKTICK, TAB, Q, TWO));
		keyToNeighborSetMap.put(TWO, getSet(ONE, Q, W, THREE));
		keyToNeighborSetMap.put(THREE, getSet(TWO, W, E, FOUR));
		keyToNeighborSetMap.put(FOUR, getSet(THREE, E, R, FIVE));
		keyToNeighborSetMap.put(FIVE, getSet(FOUR, R, T, SIX));
		keyToNeighborSetMap.put(SIX, getSet(FIVE, T, Y, SEVEN));
		keyToNeighborSetMap.put(SEVEN, getSet(SIX, Y, U, EIGHT));
		keyToNeighborSetMap.put(EIGHT, getSet(SEVEN,U,I,NINE));
		keyToNeighborSetMap.put(NINE, getSet(EIGHT, I,O,ZERO));
		keyToNeighborSetMap.put(ZERO, getSet(NINE, O,P,MINUS));
		keyToNeighborSetMap.put(MINUS, getSet(ZERO, P, LEFT_BRACKET, EQUALS));
		keyToNeighborSetMap.put(EQUALS, getSet(MINUS, LEFT_BRACKET, RIGHT_BRACKET));
		
				
		return keyToNeighborSetMap;
	}
//	
	//Helper method to help organize key-values into a set
	//useful for creating a map
	private static Set<Key> getSet(Key... keys)
	{
		return new HashSet<Key>(Arrays.asList(keys));
		
	}//end of getSet().
}//end of AppleNumericMB110LLKeyboardMetricsImpl_Luong().
