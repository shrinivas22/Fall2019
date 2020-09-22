package Examples;
public class ConwayLife {

	public static int Adder(int[][] cells, int i, int j) {
		int sum = 0;
		for (int k = i - 1; k < i + 2; k++) {
			for (int l = j - 1; l < j + 2; l++) {
				if (k >= 0 && l >= 0 && k<cells.length && l<cells[0].length) {
					if (!(k == i && l == j)) {
						sum = sum + cells[k][l];
					}
				}
			}
		}
		return sum;
	}


	public static int[][] cropper(int[][] cells) {
		int rowsumCounter1 = 0;
		int rowsumCounter2 = 0;
		int colsumCounter1 = 0;
		int colsumCounter2 = 0;
		int temp1 = 0;
		int temp2 = 0;
		for (int i = 0; i < cells.length; i++) {
			if (temp1 > 0) {
				break;
			}
			int rowsum = 0;
			for (int j = 0; j < cells[0].length; j++) {
				rowsum += cells[i][j];
			}
			if (rowsum == 0) {
				rowsumCounter1++;
			}
			temp1=rowsum;
		}
		for (int j =0; j < cells[0].length; j++) {
			if (temp2 > 0) {
				break;
			}
			int colsum = 0;
			for (int i = 0; i < cells.length; i++) {
				colsum += cells[i][j];
			}
			if (colsum == 0) {
				colsumCounter1++;
			}
			temp2=colsum;
		}
		temp1=0;temp2=0;
		for (int i =  cells.length-1 ; i>0; i--) {
			if (temp1 > 0) {
				break;
			}
			int rowsum = 0;
			for (int j = 0; j < cells[0].length; j++) {
				rowsum += cells[i][j];
			}
			if (rowsum == 0) {
				rowsumCounter2++;
			}
			temp1=rowsum;
		}
		for (int j = cells[0].length-1; j >0; j--) {
			if (temp2 > 0) {
				break;
			}
			int colsum = 0;
			for (int i = 0; i < cells.length; i++) {
				colsum += cells[i][j];
			}
			if (colsum == 0) {
				colsumCounter2++;
			}
			temp2=colsum;
		}
		//System.out.println("rw1,col1,rw2,col2:" + rowsumCounter1+","+colsumCounter1+","+rowsumCounter2+","+colsumCounter2);
		int[][] temp = new int[cells.length - rowsumCounter1 -rowsumCounter2][cells[0].length -colsumCounter1-colsumCounter2];
		for (int i = 0; i < temp.length; i++) {
			for (int j = 0; j < temp[i].length; j++) {
				temp[i][j] = cells[i+rowsumCounter1][j+colsumCounter1];
			}
		}
		return temp;
	}

public static int[][] padArray(int[][] arr, int padWith, int numOfPads) {
    int[][] temp = new int[arr.length + numOfPads*2][arr[0].length + numOfPads*2];
    for (int i = 0; i < temp.length; i++) {
        for (int j = 0; j < temp[i].length; j++) {
            temp[i][j] = padWith;
        }
    }
    for (int i = 0; i < arr.length; i++) {
        for (int j = 0; j < arr[i].length; j++) {
            temp[i+numOfPads][j+numOfPads] = arr[i][j];
        }
    }
    return temp;
}

public static int [][] getGeneration(int [][] cells, int generations)
{
//System.out.println("ORIGINAL:"+LifeDebug.htmlize(cells));
int [][] temp=  padArray(cells, 0, generations);
//System.out.println(LifeDebug.htmlize(temp));
System.out.println("length"+ temp.length+"::"+temp[0].length);
int [][] out= modifier(temp, generations);
//System.out.println(LifeDebug.htmlize(out));
return cropper(out);
}

	public static int[][] modifier(int[][] cells, int generations) {

		System.out.println("gen:" + generations);

		if (generations == 0) {
			return cells;
			// return new int[0][];
		} else {
			modifier(cells, generations - 1);
			int[][] cellsCopy = new int[cells.length][];
			for (int i = 0; i < cells.length; i++)
				cellsCopy[i] = cells[i].clone();

			for (int i = 0; i < cells.length; i++) {
				for (int j = 0; j < cells[0].length; j++) {
					int value = Adder(cellsCopy, i, j);
					if (cells[i][j] != 0) {
						if (value < 2 || value > 3) {
							cells[i][j] = 0;
						}
					} else {
						if (value == 3) {
							cells[i][j] = 1;
						}
					}
				}
			}
		}
		// your code goes here
		return cells;
	}
 	
}

