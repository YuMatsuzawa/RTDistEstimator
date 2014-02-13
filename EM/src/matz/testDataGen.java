package matz;

import java.io.*;
import java.util.Date;
import java.util.Random;

public class testDataGen {
	//EM�A���S���Y���̃e�X�g�̂��߂ɁA�Ȃ�ƂȂ��������K���z�ɏ]���Ă���悤�Ȋ����̃e�X�g�f�[�^���o�͂���B
	//Java�̐��K����(����0,�W���΍�1)���g���A������e���z���ƓK���Ȕ{��scale[k]�Ɋg��Atrans[k]���s�ړ�����B
		
	public static int cluster (int clusterMax) {//�Ă΂�邽�тɁA�ݒ肳�ꂽ���̃N���X�^�̂����A�ǂꂩ�ЂƂ̃N���X�^���w�肵�ĕԂ�
		double rnd = Math.random();
		double minval = Double.MAX_VALUE;
		int index = 0;
		for (int k=0; k<clusterMax; k++) {
			double dist = Math.abs(rnd - (double)k/(double)clusterMax);
			if (dist < minval) {
				minval = dist;
				index = k;
			}
		}
		return index;
	}
	
	public static void main (String args[]) {
		//�o�̓t�@�C���͉��s�ŋ�؂���1�����x�N�g���ŁA1�s�ɒl(���Ԃ�z��)�̂ݏ������B
		//�t�@�C���͓K���Ȑ�����悤�ɂ���B
		try {
			int fileNum = 1; //�o�̓t�@�C���̐�
			int clusterMax = 3; //�N���X�^�̐�(�����z�̕�̐�)
			if (args.length >= 1) { //����������Ȃ�o�̓t�@�C���̐��ƃN���X�^�̐��ɂ��ꂼ��Ή�������
				fileNum = Integer.parseInt(args[0]);
				if (args.length >= 2) {
					clusterMax = Integer.parseInt(args[1]);
				}
			}
			int size = 20000; 
			
			//���s���Ƀf�B���N�g�����Ȃ���΍��
			Date date = new Date();
			String dirNum = Long.toString(date.getTime());
			String path = "data" + dirNum.substring(5);
			File outPath = new File(path);
			if (!outPath.isDirectory()) {
				outPath.mkdir();
			}
			
			int itr = 1;
			int denom = 0;
			while(itr<=clusterMax){ //�d�݌���p�̕��������
				denom += itr;
				itr++;
			}
			
			for (int f = 0; f < fileNum; f++) {
				
				//double[] p = new double[clusterMax]; //�d��(�����m��)
				double[] phi = new double[clusterMax];
				double[] trans = new double[clusterMax]; //���s�ړ���
				double[] mean = new double[clusterMax];
				double[] scale = new double[clusterMax]; //�g�嗦
				double[] var = new double[clusterMax]; //���U
				int[] dataNum = new int[clusterMax];
				
				double transMin = 0.0; //���ό���p�̉����l
				for (int k=0; k<clusterMax; k++) { //���z�̃p�����[�^������
					//p[k] = (clusterMax - k)/denom;
					trans[k] = transMin + 5000*(k+1)/denom + Math.random()*1000 - 500;
					transMin = trans[k];
					scale[k] = Math.random()*2500;
					phi[k] = 0.0;
					var[k] = 0.0;
				}
				
				int[] counter = new int[100];
				for (int c=0; c<100; c++) {
					counter[c] = 0;
				}
				
				//�o��
				File outf = new File(path + String.format("/%04d.csv", f));
				OutputStreamWriter osw = new OutputStreamWriter(new FileOutputStream(outf));
				BufferedWriter bw = new BufferedWriter(osw);
				
				double[] data = new double[size];
				int[] cluster = new int[size];
				Random rnd = new Random();
				for (int i=0; i<size; i++) {//�f�[�^����
					int index = cluster(clusterMax);
					cluster[i] = index;
					dataNum[index]++;
					phi[index]++;
					//rnd.setSeed(date.getTime());
					while (true){
						double raw = rnd.nextGaussian();
						data[i] = raw*scale[index] + trans[index];
						if (data[i] >= 0 && data[i] < 10000) {
							break;
						}
					}
					int cIndex = (int)data[i]/100;
					counter[cIndex]++;
					
					mean[index] += data[i];
				}
				
				for (int k=0; k<clusterMax; k++) {
					phi[k] /= size;
					mean[k] /= (double)dataNum[k];
				}
				
				for (int i=0; i<size; i++){
					int index = cluster[i];
					double gap = data[i] - mean[index];
					var[index] += gap*gap;
				}
				
				for (int k=0; k<clusterMax; k++) {
					var[k] /= (double)dataNum[k];
				}
				
				//�f�[�^�o��
				for (int i=0; i<size; i++) {
					bw.write(String.format("%d,%f", cluster[i], data[i]));
					bw.newLine();
				}
				bw.newLine();
				
				//�x�����z�\�o��
				for (int c=0; c<100; c++) {
					bw.write(String.format("%d,%d", c*100, counter[c]));
					bw.newLine();
				}
				bw.newLine();
				
				//�p�����[�^�o��
				for(int k=0; k<clusterMax; k++) {
					bw.write(String.format("#Cluster No.%d", k));
					bw.newLine();
					bw.write(String.format("weight: %f", phi[k]));
					bw.newLine();
					bw.write(String.format("mean: %f", mean[k]));
					bw.newLine();
					bw.write(String.format("variance: %f", var[k]));
					bw.newLine();
					bw.newLine();
				}
				
				bw.close();
				osw.close();
			}
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
}