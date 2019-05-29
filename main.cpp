#include <stdio.h>
#include <iostream>
#include <time.h>
#include <direct.h>

#include <pcl/visualization/cloud_viewer.h>
#include <pcl/io/io.h>
#include <pcl/io/pcd_io.h>
#include <pcl/io/ply_io.h>

#include <vector>
#include <opencv2/opencv.hpp>
#include <Eigen/Core>

std::string DIR = "gene/texture";
bool RGBSET_MODE = false;
bool INTERRUPT_MODE = false;
bool OCLUSION_OUTPUT = false;
std::string IN_FILENAME = "mesh";
std::string OUT_FILENAME = "mesh3d";
double CAMERAOFST_V = 0.0;
double CAMERAOFST_PAN = 0.0;
double CAMERAOFST_TILT = 0.0;
double CAMERAOFST_ROLL = 0.0;
double UVMAP_OFST_U = 0;
double UVMAP_OFST_V = 0;
double UVMAP_RATIO = 1.0;
double UVMAP_TRIM = 1.0;
double GROUND_ANGLE_TH = 40.0;
double INCIDENCE_ANGLE_TH = 85.0;
double DISTANCE_TH = 15.0;
double OC_ANGLE_TH = 3.0;

std::vector<std::string> getImageName(std::string dir, std::string ini) {
	HANDLE hFind;
	WIN32_FIND_DATA win32fd;
	std::vector<std::string> file_names;
	std::string search_name = dir + "/" + ini + "*.jpg";
	hFind = FindFirstFile(search_name.c_str(), &win32fd);
	do {
		if (win32fd.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY) {}
		else {
			file_names.push_back(win32fd.cFileName);
		}
	} while (FindNextFile(hFind, &win32fd));
	FindClose(hFind);
	return file_names;
}

int setParam() {
	// ���� �p�����[�^�Ǎ���
	// --------------------------------------------------------------------------------------------
	std::ifstream fin("parameter.txt");
	if (!fin) return 1;
	std::string line;
	while (std::getline(fin, line))
	{
		if (line[0] == '#') continue;
		if (line.find('=') == std::string::npos) continue;
		std::stringstream ss(line);
		std::string name;
		ss >> name;
		ss.ignore(line.size(), '=');
		if (name == "TextureMapping.DIR") { ss >> DIR; }
		else if (name == "TextureMapping.RGBSET_MODE") { ss >> RGBSET_MODE; }
		else if (name == "TextureMapping.INTERRUPT_MODE") { ss >> INTERRUPT_MODE; }
		else if (name == "TextureMapping.OCLUSION_OUTPUT") { ss >> OCLUSION_OUTPUT; }
		else if (name == "TextureMapping.IN_FILENAME") { ss >> IN_FILENAME; }
		else if (name == "TextureMapping.OUT_FILENAME") { ss >> OUT_FILENAME; }
		else if (name == "TextureMapping.CAMERAOFST_V") { ss >> CAMERAOFST_V; }
		else if (name == "TextureMapping.CAMERAOFST_PAN") { ss >> CAMERAOFST_PAN; }
		else if (name == "TextureMapping.CAMERAOFST_TILT") { ss >> CAMERAOFST_TILT; }
		else if (name == "TextureMapping.CAMERAOFST_ROLL") { ss >> CAMERAOFST_ROLL; }
		else if (name == "TextureMapping.UVMAP_OFST_U") { ss >> UVMAP_OFST_U; }
		else if (name == "TextureMapping.UVMAP_OFST_V") { ss >> UVMAP_OFST_V; }
		else if (name == "TextureMapping.UVMAP_RATIO") { ss >> UVMAP_RATIO; }
		else if (name == "TextureMapping.UVMAP_TRIM") { ss >> UVMAP_TRIM; }
		else if (name == "TextureMapping.GROUND_ANGLE_TH") { ss >> GROUND_ANGLE_TH; }
		else if (name == "TextureMapping.INCIDENCE_ANGLE_TH") { ss >> INCIDENCE_ANGLE_TH; }
		else if (name == "TextureMapping.DISTANCE_TH") { ss >> DISTANCE_TH; }
		else if (name == "TextureMapping.OC_ANGLE_TH") { ss >> OC_ANGLE_TH; }
	}

	std::cout << "TextureMapping.DIR_______________: " << DIR << std::endl;
	std::cout << "TextureMapping.RGBSET_MODE_______: " << RGBSET_MODE << std::endl;
	std::cout << "TextureMapping.INTERRUPT_MODE____: " << INTERRUPT_MODE << std::endl;
	std::cout << "TextureMapping.OCLUSION_OUTPUT___: " << OCLUSION_OUTPUT << std::endl;
	std::cout << "TextureMapping.IN_FILENAME_______: " << IN_FILENAME << std::endl;
	std::cout << "TextureMapping.OUT_FILENAME______: " << OUT_FILENAME << std::endl;
	std::cout << "TextureMapping.CAMERAOFST_V______: " << CAMERAOFST_V << std::endl;
	std::cout << "TextureMapping.CAMERAOFST_PAN____: " << CAMERAOFST_PAN << std::endl;
	std::cout << "TextureMapping.CAMERAOFST_TILT___: " << CAMERAOFST_TILT << std::endl;
	std::cout << "TextureMapping.CAMERAOFST_ROLL___: " << CAMERAOFST_ROLL << std::endl;
	std::cout << "TextureMapping.UVMAP_OFST_U______: " << UVMAP_OFST_U << std::endl;
	std::cout << "TextureMapping.UVMAP_OFST_V______: " << UVMAP_OFST_V << std::endl;
	std::cout << "TextureMapping.UVMAP_RATIO_______: " << UVMAP_RATIO << std::endl;
	std::cout << "TextureMapping.UVMAP_TRIM________: " << UVMAP_TRIM << std::endl;
	std::cout << "TextureMapping.GROUND_ANGLE_TH___: " << GROUND_ANGLE_TH << std::endl;
	std::cout << "TextureMapping.INCIDENCE_ANGLE_TH: " << INCIDENCE_ANGLE_TH << std::endl;
	std::cout << "TextureMapping.DISTANCE_TH_______: " << DISTANCE_TH << std::endl;
	std::cout << "TextureMapping.OC_ANGLE_TH_______: " << OC_ANGLE_TH << std::endl;

	return 0;
}

int main()
{
	//////////////01234567890123456789012345678901234567890123456789
	std::cout << "--------------------------------------------------" << std::endl;
	std::cout << "- Start Texture Mapping                           " << std::endl;
	std::cout << "--------------------------------------------------" << std::endl;

	if (setParam()) {
		std::cout << "- Pamameter Read Error!!! " << std::endl;
	};

	// �t�@�C�����ݒ�
	std::string in_file_ply = IN_FILENAME + ".ply";
	std::string out_file_ply = OUT_FILENAME + ".ply";		//�o��PLY�t�@�C����
	std::string out_file_obj = OUT_FILENAME + ".obj";		//�o��OBJ�t�@�C����
	std::string out_file_mtl = OUT_FILENAME + ".mtl";		//�o��MTL�t�@�C����
	std::string out_file_texture = OUT_FILENAME + ".jpg";	//�o�̓e�N�X�`����

	// ���� �摜�t�@�C���Ǎ���
	// --------------------------------------------------------------------------------------------
	std::vector<std::string> file_names_a = getImageName(DIR, "data_a");
	std::vector<std::string> file_names_b = getImageName(DIR, "data_b");

	//for (auto f : file_names) {std::cout << f << std::endl;}
	int tx_num = file_names_a.size();						//�e�N�X�`���t�@�C����(�v����)
	std::sort(file_names_a.begin(), file_names_a.end());	//�t�@�C������(�擾��)�Ƀ\�[�g
	std::sort(file_names_b.begin(), file_names_b.end());	//�t�@�C������(�擾��)�Ƀ\�[�g
	for (int i = 0; i < tx_num; i++) {
		std::cout << "File_a[" << i << "] : " << file_names_a[i] << std::endl;
		std::cout << "File_b[" << i << "] : " << file_names_b[i] << std::endl;
	}
	// ���̓e�N�X�`���Ǎ��݁C�o�̓e�N�X�`���̈�ݒ�
	std::vector<cv::Mat> src_img_a, src_img_b;

	for (int i = 0; i < tx_num; i++) {
		cv::Mat src_a = cv::imread(DIR + "/" + file_names_a[i], 1);	//�O���摜�Ǎ�
		//cv::Mat dst_a;
		cv::resize(src_a, src_a, cv::Size(), UVMAP_RATIO, UVMAP_RATIO); //�w��{���ŉ𑜓x�ύX
		src_img_a.push_back(src_a);

		cv::Mat src_b = cv::imread(DIR + "/" + file_names_b[i], 1);	//�㔼�摜�Ǎ�
		//cv::Mat dst_b;
		src_b = cv::imread(DIR + "/" + file_names_b[i], 1);	//�摜�Ǎ�
		cv::resize(src_b, src_b, cv::Size(), UVMAP_RATIO, UVMAP_RATIO); //�w��{���ŉ𑜓x�ύX
		src_img_b.push_back(src_b);
	}

	cv::Mat link_img, link_img_a, link_img_b;
	cv::vconcat(src_img_a, link_img_a);				//�O���摜�F�c�����ɘA���`
	cv::vconcat(src_img_b, link_img_b);				//�㔼�摜�F�c�����ɘA���a
	cv::hconcat(link_img_a, link_img_b, link_img);	//�`�Ƃa���������ɘA��
	std::cout << "Texture Size (x , y): " << link_img.cols << " , " << link_img.rows << std::endl;

	// �����@�v���n�_-�ϊ��s��Ǎ��݁C�t�s��v�Z ���ő�30�ӏ�
	// --------------------------------------------------------------------------------------------
	//////////////01234567890123456789012345678901234567890123456789
	std::cout << "- Read Transformation matrix ---------------------" << endl;
	Eigen::Matrix4d trans0[30], inv0[30];	//���ʌX���␳
	Eigen::Matrix4d trans1[30], inv1[30];	//�v���ʒu�ϊ��i����v���ʒu�����_�Ƃ���j
	ifstream fin0(DIR + "/trans0s.dat");
	ifstream fin1(DIR + "/trans1s.dat");
	for (int k = 0; k < tx_num; k++) {
		for (int i = 0; i < 4; i++) {
			for (int j = 0; j < 4; j++) {
				fin0 >> trans0[k](i, j);
				fin1 >> trans1[k](i, j);
			}
		}
		inv0[k] = trans0[k].inverse();
		inv1[k] = trans1[k].inverse();
	}

	// �����@PLY�t�@�C���Ǎ���
	// --------------------------------------------------------------------------------------------
	//////////////01234567890123456789012345678901234567890123456789
	std::cout << "- Read Mesh Ply Data -----------------------------" << endl;
	pcl::PLYReader reader;
	pcl::PolygonMesh in_mesh;
	reader.read(DIR + "/" + in_file_ply, in_mesh);
	pcl::PointCloud<pcl::PointXYZ> in_cloud;
	pcl::fromPCLPointCloud2(in_mesh.cloud, in_cloud);		//PLY�̓_�Q��PCD��
	//////////////012345678901234567890123456789
	std::cout << "Vertex Num____________________: " << in_mesh.cloud.width << endl;		//���_��
	std::cout << "Face Num______________________: " << in_mesh.polygons.size() << endl; //�ʐ�

	// �����@���b�V���̖@���x�N�g��
	// --------------------------------------------------------------------------------------------
	const int m_dims = 3;	//���_���� (xyz)
	std::vector<pcl::PointXYZ> pm0, pn0;
	for (int i = 0; i < in_mesh.polygons.size(); i++) {
		pcl::PointXYZ p[m_dims];
		pcl::PointXYZ pm, pn;
		for (int j = 0; j < m_dims; j++) {
			p[j] = in_cloud.points[in_mesh.polygons[i].vertices[j]];
		}
		//���b�V���̖ʒ��S ---------------------------------------------------
		pm.x = (p[0].x + p[1].x + p[2].x) / 3.0;
		pm.y = (p[0].y + p[1].y + p[2].y) / 3.0;
		pm.z = (p[0].z + p[1].z + p[2].z) / 3.0;
		pm0.push_back(pm);
		//���b�V���̖@���x�N�g���v�Z (�Q�Ӄx�N�g���̊O��) --------------------
		float p01x = p[1].x - p[0].x;
		float p01y = p[1].y - p[0].y;
		float p01z = p[1].z - p[0].z;
		float p12x = p[2].x - p[1].x;
		float p12y = p[2].y - p[1].y;
		float p12z = p[2].z - p[1].z;
		pn.x = p01y * p12z - p01z * p12y;
		pn.y = p01z * p12x - p01x * p12z;
		pn.z = p01x * p12y - p01y * p12x;
		pn0.push_back(pn);
	}
	//////////////012345678901234567890123456789
	std::cout << "Mesh Normal Vector Num________: " << pn0.size() << endl; //�ʐ�

	// ���� ���b�V���ƃe�N�X�`���̑Ή��t��
	// --------------------------------------------------------------------------------------------
	//////////////01234567890123456789012345678901234567890123456789
	std::cout << "- Start Mapping Calcuration ----------------------" << std::endl;

	// ���@�J�����ƃf�W�^�C�U�̊p�x�␳(Pan:z�CTilt:x Roll:y)
	const float	p_rad = DEG2RAD(CAMERAOFST_PAN);
	const float	t_rad = DEG2RAD(CAMERAOFST_TILT);
	const float	r_rad = DEG2RAD(CAMERAOFST_ROLL);

	//�v���_-���b�V�����S�ԃx�N�g���ƁC���b�V���@���x�N�g���̊p�x���{�p�x�ȏ�̏ꍇ�̓e�N�X�`����\��Ȃ��D
	//�v���_���猩�Ėʂ������Ȃ��C�܂��C�s�p�̏ꍇ�̓e�N�X�`����\��Ȃ��D
	const float nr_rad = DEG2RAD(INCIDENCE_ANGLE_TH);

	//�f�W�^�C�U�̉摜���������邽�߁C�v���_�ƃO�����h�ԃx�N�g���Ɩ{�p�x�͈̔͂̏��ʃ��b�V���ɂ̓��b�V����\��Ȃ��D
	const float gnd_rad = DEG2RAD(GROUND_ANGLE_TH);

	const int t_dims = 6;	//UV�}�b�v�w�莟��:uv(2)�~���_��(3)
	float uv[t_dims];
	int count = 0;
	int mesh_num = 0;
	std::vector<bool> mesh_enable;
	std::vector<pcl::PointXY> uv1, uv2, uv3;
	for (int i = 0; i < in_mesh.polygons.size(); i++) {

		pcl::PointXYZ p[m_dims];

		//���b�V�����_�擾
		for (int j = 0; j < m_dims; j++) {
			p[j] = in_cloud.points[in_mesh.polygons[i].vertices[j]];
		}

		//���b�V�����S�_�v�Z ------------------------------------------------
		float pmx = (p[0].x + p[1].x + p[2].x) / 3.0;
		float pmy = (p[0].y + p[1].y + p[2].y) / 3.0;
		float pmz = (p[0].z + p[1].z + p[2].z) / 3.0;

		//���b�V���@���x�N�g���v�Z (�Q�Ӄx�N�g���̊O��)----------------------
		float p01x = p[1].x - p[0].x;
		float p01y = p[1].y - p[0].y;
		float p01z = p[1].z - p[0].z;
		float p12x = p[2].x - p[1].x;
		float p12y = p[2].y - p[1].y;
		float p12z = p[2].z - p[1].z;
		float pnx = p01y * p12z - p01z * p12y;
		float pny = p01z * p12x - p01x * p12z;
		float pnz = p01x * p12y - p01y * p12x;

		//���b�V���ɋ߂��v���ʒu(tx_id)�̎擾�D
		int tx_id = 0;
		bool tx_enable = false;
		float lmax = DISTANCE_TH; // scan_range_max;
		for (int k = 0; k < tx_num; k++) {
			float dx = trans1[k](0, 3) - pmx;
			float dy = trans1[k](1, 3) - pmy;
			float dz = trans1[k](2, 3) - pmz;
			float lxy = sqrt(dx*dx + dy * dy);
			float direc = (pnx*dx + pny * dy + pnz * dz) / sqrt(pnx*pnx + pny * pny + pnz * pnz) / sqrt(dx*dx + dy * dy + dz * dz);
			// ���ˊp�ɂ�锻��
			if ((direc > cos(nr_rad)) && (lxy / dz > tan(gnd_rad) || dz < 0)) {
				float l = sqrt(dx*dx + dy * dy + dz * dz);

				//�����ɂ�锻�� -----------------------------------------------------
				if (l < lmax) {
					bool interrupt = false;

					// �Օ����̉e�� ��True�ɂ���ƕ��ב�-------------------------------
					if (INTERRUPT_MODE) {

						for (int m = 0; m < in_mesh.polygons.size(); m++) {

							// �v���_����Q�ƃ��b�V�����S�܂ł̋���:ml
							float mdx = trans1[k](0, 3) - pm0[m].x;
							float mdy = trans1[k](1, 3) - pm0[m].y;
							float mdz = trans1[k](2, 3) - pm0[m].z;
							float ml = sqrt(mdx*mdx + mdy * mdy + mdz * mdz);

							if (ml < l) {
								// �Q�ƃ��b�V������{���b�V�����S�܂ł̋���:dl
								float ddx = pm0[m].x - pmx;
								float ddy = pm0[m].y - pmy;
								float ddz = pm0[m].z - pmz;
								float dl = sqrt(ddx*ddx + ddy * ddy + ddz * ddz);

								if (dl < l) {
									float md_direc = (ddx*dx + ddy * dy + ddz * dz) / sqrt(ddx*ddx + ddy * ddy + ddz * ddz) / sqrt(dx*dx + dy * dy + dz * dz);
									if (md_direc > cos(DEG2RAD(OC_ANGLE_TH))) {
										interrupt = true;
										break;
									}
								}
							}
						}
					}

					if (!interrupt) { 		// �����ɂ�锻��
						lmax = l;
						tx_id = k + 1;
					}
				}
			}
		}

		// ���@�f�W�^�C�U�̌v���ʒu�Ɋ�Â��ʒu�␳
		if (tx_id != 0) {
			for (int j = 0; j < m_dims; j++) {
				Eigen::MatrixXd Pin(4, 1), Pout(4, 1);
				Pin << p[j].x, p[j].y, p[j].z, 1;
				Pout = inv0[tx_id - 1] * inv1[tx_id - 1] * Pin;
				p[j].x = Pout(0);
				p[j].y = Pout(1);
				p[j].z = Pout(2);
			}
		}

		// ���@�J�����|�f�W�^�C�U�Ԏp���␳
		for (int j = 0; j < m_dims; j++) {
			Eigen::MatrixXd Pin(3, 1), Rp(3, 3), Rr(3, 3), Rt(3, 3), Pout(3, 1);
			Pin << p[j].x, p[j].y, p[j].z - CAMERAOFST_V;
			Rp << cos(p_rad), -sin(p_rad), 0, sin(p_rad), cos(p_rad), 0, 0, 0, 1;

			Rr << cos(r_rad), 0, sin(r_rad), 0, 1, 0, -sin(r_rad), 0, cos(r_rad);

			Rt << 1, 0, 0, 0, cos(t_rad), -sin(t_rad), 0, sin(t_rad), cos(t_rad);
			Pout = Rt * Rr * Rp * Pin;
			p[j].x = Pout(0);
			p[j].y = Pout(1);
			p[j].z = Pout(2);
		}

		// ���@���b�V�����_�ɑ�������UV�}�b�v�ʒu�Z�o,�e�N�X�`���̗L����������
		if (tx_id != 0) {
			mesh_enable.push_back(true);
			mesh_num++;
			for (int j = 0; j < m_dims; j++) {

				// -----�܂�f�W�̉摜�ƈقȂ鏈�� ------------------------------------
				double l = sqrt(p[j].x * p[j].x + p[j].y * p[j].y + p[j].z * p[j].z);
				double l_yz = sqrt(p[j].y * p[j].y + p[j].z * p[j].z);
				double ang_l = asin(l_yz / l) * UVMAP_TRIM;
				
				double ang_yz = atan2(p[j].z, p[j].y);
				double ang_v = ang_l * sin(ang_yz);
				double ang_h = ang_l * cos(ang_yz);
				//double ang_v = ang_l * sin(p[j].y / l_yz);
				//double ang_h = ang_l * cos(p[j].z / l_yz);

				uv[j * 2 + 1] = (ang_v / M_PI + tx_num - tx_id + 0.5f + UVMAP_OFST_V) / tx_num;
				if (uv[j * 2 + 1] < 0) { uv[j * 2 + 1] = 0; }
				else if (uv[j * 2 + 1] > 1) { uv[j * 2 + 1] = 1; }

				if (p[0].x > 0) {
					ang_h = M_PI / 2.0f - ang_h;
				}
				else {
					ang_h = 3.0 * M_PI / 2.0f + ang_h;
				}
				uv[j * 2] = ang_h / M_PI / 2.0f + UVMAP_OFST_U;
				if (uv[j * 2] < 0) { uv[j * 2] = 0; }
				else if (uv[j * 2] > 1) { uv[j * 2] = 1; }
				// --------------------------------------------------------------------

			}
			if (p[0].y < 0) {
				if (p[0].x < 0 && (p[1].x > 0 || p[2].x > 0)) { uv[0] = uv[0]; }
				if (p[1].x < 0 && (p[2].x > 0 || p[0].x > 0)) { uv[2] = uv[2]; }
				if (p[2].x < 0 && (p[0].x > 0 || p[1].x > 0)) { uv[4] = uv[4]; }
			}
		}
		else {	// ���b�V���ɑ�������UV�}�b�v�������ꍇ�C�摜-�����̐F�����蓖�Ă�D
			mesh_enable.push_back(false);
			uv[0] = 0; uv[1] = 0;
			uv[2] = 0; uv[3] = 0.02/tx_num; //uv[2] = 0; uv[3] = 10.0 / link_img.size().height;
			uv[4] = 0.02; uv[5] = 0; //uv[4] = 10.0 / link_img.size().width; uv[5] = 0;
		}

		pcl::PointXY p1, p2, p3;
		p1.x = uv[0];
		p1.y = uv[1];
		p2.x = uv[2];
		p2.y = uv[3];
		p3.x = uv[4];
		p3.y = uv[5];
		uv1.push_back(p1);
		uv2.push_back(p2);
		uv3.push_back(p3);

		// ���@���k�Ńe�N�X�`���FUV�}�b�v �؎�C�\�t
		float u_min = uv[0], u_max = uv[0];
		float v_min = uv[1], v_max = uv[1];
		for (int j = 1; j < m_dims; j++) {
			if (uv[j * 2] < uv[j * 2 - 2] && uv[j * 2] < u_min) { u_min = uv[j * 2]; }
			if (uv[j * 2] > uv[j * 2 - 2] && uv[j * 2] > u_max) { u_max = uv[j * 2]; }
			if (uv[j * 2 + 1] < uv[j * 2 - 1] && uv[j * 2 + 1] < v_min) { v_min = uv[j * 2 + 1]; }
			if (uv[j * 2 + 1] > uv[j * 2 - 1] && uv[j * 2 + 1] > v_max) { v_max = uv[j * 2 + 1]; }
		}
		int img_x = int(u_min * link_img.cols);
		int img_dx = int(u_max * link_img.cols) - img_x;
		int img_y = int((1.0 - v_max) * link_img.rows);
		int img_dy = int((1.0 - v_min) * link_img.rows) - img_y;
		if (img_x + img_dx < link_img.cols && img_y + img_dy < link_img.rows && tx_id != 0) {
			cv::Mat cut_img(link_img, cv::Rect(img_x, img_y, img_dx, img_dy));
			cv::Mat roi = link_img(cv::Rect(img_x, img_y, img_dx, img_dy));
			cut_img.copyTo(roi);
		}

		// ���@�i���J�E���^�\��(1000���b�V���P��)
		if (count == 10000) {
			std::cout << "Mesh Calc Progress : " << i << "/" << in_mesh.polygons.size() << endl;
			count = 1;
		}
		else { count++; }
	}

	// �����@�o��
	// --------------------------------------------------------------------------------------------

	// ���@UV�}�b�v�o��(.jpg)
	cv::imwrite(DIR + "/" + out_file_texture, link_img);
	// --------------------------------------------------------------------------------------------

	// �� �_�Q(data_n.ply)�ւ�RGB�ǉ�
	if (RGBSET_MODE) {
		//////////////01234567890123456789012345678901234567890123456789
		std::cout << "- Start RGB Set ----------------------------------" << std::endl;

		//*
		// �f�[�^�ǂݍ���
		//////////////01234567890123456789012345678901234567890123456789
		std::cout << "- Ply Data -----------------------------" << endl;
		pcl::PointCloud<pcl::PointNormal> ::Ptr in_cloud_n(new pcl::PointCloud<pcl::PointNormal>);
		reader.read(DIR + "/data_n.ply", *in_cloud_n);
		pcl::PointCloud<pcl::PointXYZRGB>::Ptr out_cloud_rgb(new pcl::PointCloud<pcl::PointXYZRGB>);

		count = 0;
		for (int i = 0; i < in_cloud_n->size(); i++) {
			pcl::PointXYZRGB p_rgb;
			pcl::PointXYZ p;
			p_rgb.x = in_cloud_n->points[i].x;
			p_rgb.y = in_cloud_n->points[i].y;
			p_rgb.z = in_cloud_n->points[i].z;

			//�摜�ʒu(tx_id)�̎擾�D
			int tx_id = 0;
			float lmax = DISTANCE_TH; // scan_range_max;
			for (int k = 0; k < tx_num; k++) {
				float dx = trans1[k](0, 3) - p_rgb.x;
				float dy = trans1[k](1, 3) - p_rgb.y;
				float dz = trans1[k](2, 3) - p_rgb.z;
				//			float lxy = sqrt(dx*dx + dy*dy);
				float l = sqrt(dx*dx + dy * dy + dz * dz);

				//�����ɂ�锻�� -----------------------------------------------------

				if (l < lmax) {
					lmax = l;
					tx_id = k + 1;
				}
			}

			// �f�W�^�C�U�̌v���ʒu�Ɋ�Â��ʒu�␳
			if (tx_id != 0) {
				Eigen::MatrixXd Pin(4, 1), Pout(4, 1);
				Pin << p_rgb.x, p_rgb.y, p_rgb.z, 1;
				Pout = inv0[tx_id - 1] * inv1[tx_id - 1] * Pin;
				p.x = Pout(0);
				p.y = Pout(1);
				p.z = Pout(2);
			}

			// �J�����|�f�W�^�C�U�Ԏp���␳
			Eigen::MatrixXd Pin(3, 1), Rp(3, 3), Rr(3, 3), Rt(3, 3), Pout(3, 1);
			Pin << p.x, p.y, p.z - CAMERAOFST_V;
			Rp << cos(p_rad), -sin(p_rad), 0, sin(p_rad), cos(p_rad), 0, 0, 0, 1;
			Rr << cos(r_rad), 0, sin(r_rad), 0, 1, 0, -sin(r_rad), 0, cos(r_rad);
			Rt << 1, 0, 0, 0, cos(t_rad), -sin(t_rad), 0, sin(t_rad), cos(t_rad);
			Pout = Rt * Rr * Rp * Pin;
			p.x = Pout(0);
			p.y = Pout(1);
			p.z = Pout(2);

			// ���b�V�����_�ɑ�������UV�}�b�v�ʒu�Z�o,�e�N�X�`���̗L����������
			//float ratio_x = (M_PI + atan2(p.x, p.y)) / (2 * M_PI);
			float ratio_x = asin(p.z / sqrt(p.z * p.z + p.x * p.x)) / M_PI;
			int img_x = int(link_img.cols * ratio_x);
			// float ratio_y = atan(p.z / sqrt(p.x * p.x + p.y * p.y)) / M_PI;
			float ratio_y = asin(p.y / sqrt(p.y * p.y + p.x * p.x)) / M_PI;
			int img_y = int(link_img.rows * (tx_id - ratio_y - 0.5f) / tx_num);

			cv::Mat3b dotImg = link_img;
			cv::Vec3b bgr = dotImg(cv::Point(img_x, img_y));
			p_rgb.b = bgr[0];
			p_rgb.g = bgr[1];
			p_rgb.r = bgr[2];
			out_cloud_rgb->push_back(p_rgb);

			// �i���J�E���^�\��(1000���b�V���P��)
			if (count == 10000) {
				std::cout << "Set RGB Calc Progress : " << i << "/" << in_cloud_n->size() << " x:" << img_x << "y:" << img_y << endl;
				count = 1;
			}
			else { count++; }
		}
		pcl::PLYWriter writer;
		writer.write(DIR + "/data_c.ply", *out_cloud_rgb, true);
		std::cout << "Saved [data_c.ply]" << std::endl;
		pcl::io::savePCDFile(DIR + "/data_c.pcd", *out_cloud_rgb, true);
		std::cout << "Saved [data_c.pcd]" << std::endl;
	}

	// ���@PLY�f�[�^�o��(.ply)
	// --------------------------------------------------------------------------------------------
	std::string fname;
	fname = DIR + "/" + out_file_ply;
	
	std::ofstream ofp(fname);
	ofp << "ply" << std::endl;
	ofp << "format ascii 1.0" << std::endl;
	ofp << "comment PCL generated" << std::endl;
	ofp << "comment TextureFile " << out_file_texture << std::endl;
	ofp << "element vertex " << in_mesh.cloud.width << std::endl;
	ofp << "property float x" << std::endl;
	ofp << "property float y" << std::endl;
	ofp << "property float z" << std::endl;
	if (OCLUSION_OUTPUT) {
		ofp << "element face " << in_mesh.polygons.size() << std::endl;
	}
	else {
		ofp << "element face " << mesh_num << std::endl;
	}
	ofp << "property list uchar int vertex_indices" << std::endl;
	ofp << "property list uchar float texcoord" << std::endl;
	ofp << "end_header" << std::endl;

	// --- ���_�o�� ---
	for (int i = 0; i < in_mesh.cloud.width; i++) {
		pcl::PointXYZ p = in_cloud.points[i];
		ofp << p.x << " " <<p.y << " " << p.z << std::endl;
	}

	// --- UV�}�b�v���� ---
	for (int i = 0; i < in_mesh.polygons.size(); i++) {
		if (mesh_enable[i] || OCLUSION_OUTPUT) {
			ofp << m_dims << " " << std::endl;
			for (int j = 0; j < m_dims; j++) {
				ofp << in_mesh.polygons[i].vertices[j] << " ";
			}
			ofp << t_dims << " " << uv1[i].x << " " << uv1[i].y << " " << uv2[i].x << " " << uv2[i].y << " " << uv3[i].x << " " << uv3[i].y << std::endl;
		}
	}

	/*
	FILE *fp_ply = fopen(fname.c_str(), "w");

	// --- �w�b�_�o�� ---
	fprintf(fp_ply, "ply\n");
	fprintf(fp_ply, "format ascii 1.0\n");
	fprintf(fp_ply, "comment PCL generated\n");	//fprintf(fp, "comment VCGLIB generated\n");
	fprintf(fp_ply, "comment TextureFile %s\n", out_file_texture.c_str());
	fprintf(fp_ply, "element vertex %d\n", in_mesh.cloud.width);//���_��
	fprintf(fp_ply, "property float x\n");
	fprintf(fp_ply, "property float y\n");
	fprintf(fp_ply, "property float z\n");
	if (OCLUSION_OUTPUT) {
		fprintf(fp_ply, "element face %d\n", in_mesh.polygons.size());//���b�V����
	}
	else {
		fprintf(fp_ply, "element face %d\n", mesh_num);//���b�V����
	}
	fprintf(fp_ply, "property list uchar int vertex_indices\n");
	fprintf(fp_ply, "property list uchar float texcoord\n");
	fprintf(fp_ply, "end_header\n");

	// --- ���_�o�� ---
	for (int i = 0; i < in_mesh.cloud.width; i++) {
		pcl::PointXYZ p = in_cloud.points[i];
		fprintf(fp_ply, "%f %f %f\n", p.x, p.y, p.z);
	}

	// --- UV�}�b�v���� ---
	for (int i = 0; i < in_mesh.polygons.size(); i++) {
		if (mesh_enable[i] || OCLUSION_OUTPUT) {
			fprintf(fp_ply, "%d ", m_dims);
			for (int j = 0; j < m_dims; j++) {
				fprintf(fp_ply, "%d ", in_mesh.polygons[i].vertices[j]);
			}
			fprintf(fp_ply, "%d %f %f %f %f %f %f \n", t_dims, uv1[i].x, uv1[i].y, uv2[i].x, uv2[i].y, uv3[i].x, uv3[i].y);
		}
	}
	fclose(fp_ply);
	*/
	//////////////01234567890123456789012345678901234567890123456789
	std::cout << "- Output PLY -------------------------------------" << endl;

	// ���f�[�^�o��(MTL)
	// --------------------------------------------------------------------------------------------
	fname = DIR + "/" + out_file_mtl;

	std::ofstream ofm(fname);
	ofm << "# Generated by 3D-Digitizer" << std::endl;
	ofm << "newmtl default" << std::endl;
	ofm << "#Ka 0.17 0.17 0.17" << std::endl;
	ofm << "#Kd 1 1 1" << std::endl;
	ofm << "Ks 0 0 0" << std::endl;
	ofm << "map_Kd " << out_file_texture << std::endl;
	/*
	FILE *fp_mtl = fopen(fname.c_str(), "w");
	fprintf(fp_mtl, "# Generated by 3D-Digitizer\n");
	fprintf(fp_mtl, "newmtl default\n");
	fprintf(fp_mtl, "#Ka 0.17 0.17 0.17\n");
	fprintf(fp_mtl, "#Kd 1 1 1\n");
	fprintf(fp_mtl, "Ks 0 0 0\n");
	fprintf(fp_mtl, "map_Kd %s\n", out_file_texture.c_str());
	fclose(fp_mtl);
	*/
	//////////////01234567890123456789012345678901234567890123456789
	std::cout << "- Output MTL -------------------------------------" << endl;

	// ���@�f�[�^�o��(OBJ)
	// --------------------------------------------------------------------------------------------
	fname = DIR + "/" + out_file_obj;

	std::ofstream ofo(fname);
	ofo << "# Generated by 3D-Digitizer" << std::endl;
	ofo << "mtllib " << out_file_mtl << std::endl;
	ofo << "o Mesh" << std::endl;
	// ���_�o��
	for (int i = 0; i < in_mesh.cloud.width; i++) {
		pcl::PointXYZ p = in_cloud.points[i];
		ofo << "v " << p.x << " " << p.y << " " << p.z << std::endl;
	}
	// UV�}�b�v����
	for (int i = 0; i < in_mesh.polygons.size(); i++) {
		if (mesh_enable[i] || OCLUSION_OUTPUT) {
			ofo << "vt " << uv1[i].x << " " << uv1[i].y << std::endl;
			ofo << "vt " << uv2[i].x << " " << uv2[i].y << std::endl;
			ofo << "vt " << uv3[i].x << " " << uv3[i].y << std::endl;
		}
	}
	// ���b�V������
	ofo << "usemtl default" << std::endl;
	count = 0;
	for (int i = 0; i < in_mesh.polygons.size(); i++) {
		if (mesh_enable[i] || OCLUSION_OUTPUT) {
			ofo << "f" ;
			for (int j = 0; j < m_dims; j++) {
				ofo << " " << in_mesh.polygons[i].vertices[j] + 1 << "/" << count*m_dims + j + 1 << "/" << in_mesh.polygons[i].vertices[j] + 1;
			}
			count++;
			ofo << std::endl;
		}
	}
	if (OCLUSION_OUTPUT) {
		ofo << "# " << in_mesh.polygons.size() << " faces" << std::endl;
	}
	else {
		ofo << "# " << mesh_num << " faces" << std::endl;
	}

	/*
	FILE *fp_obj = fopen(fname.c_str(), "w");
	fprintf(fp_obj, "# Generated by 3D-Digitizer\n");
	fprintf(fp_obj, "mtllib %s\n", out_file_mtl.c_str());
	fprintf(fp_obj, "o Mesh\n");

	// ���_�o��
	for (int i = 0; i < in_mesh.cloud.width; i++) {
		pcl::PointXYZ p = in_cloud.points[i];
		fprintf(fp_obj, "v %f %f %f\n", p.x, p.y, p.z);
	}

	// UV�}�b�v����
	for (int i = 0; i < in_mesh.polygons.size(); i++) {
		if (mesh_enable[i] || OCLUSION_OUTPUT) {
			fprintf(fp_obj, "vt %f %f\n", uv1[i].x, uv1[i].y);
			fprintf(fp_obj, "vt %f %f\n", uv2[i].x, uv2[i].y);
			fprintf(fp_obj, "vt %f %f\n", uv3[i].x, uv3[i].y);
		}
	}

	// ���b�V������
	fprintf(fp_obj, "usemtl default\n");
	count = 0;
	for (int i = 0; i < in_mesh.polygons.size(); i++) {
		if (mesh_enable[i] || OCLUSION_OUTPUT) {
			fprintf(fp_obj, "f");
			for (int j = 0; j < m_dims; j++) {
				fprintf(fp_obj, " %d/%d/%d", in_mesh.polygons[i].vertices[j] + 1, count*m_dims + j + 1, in_mesh.polygons[i].vertices[j] + 1);
			}
			count++;
			fprintf(fp_obj, "\n");
		}
	}
	if (OCLUSION_OUTPUT) {
		fprintf(fp_obj, "# %d faces\n", in_mesh.polygons.size());
	}
	else {
		fprintf(fp_obj, "# %d faces\n", mesh_num);
	}
	fclose(fp_obj);
	*/
	//////////////01234567890123456789012345678901234567890123456789
	std::cout << "- Output OBJ -------------------------------------" << endl;

	// ���� �I��
	//////////////01234567890123456789012345678901234567890123456789
	std::cout << "--------------------------------------------------" << endl;
	std::cout << "- Finish Tecture Mapping                          " << endl;
	std::cout << "--------------------------------------------------" << endl;

	return 0;
}