import sys
import math
import matplotlib
import pandas as pd
import matplotlib.pyplot as plt
from PyQt5.QtCore import QTimer, QMutex, Qt
from PyQt5.QtWidgets import QApplication, QMainWindow, QPushButton, QLabel, QSpinBox, QFileDialog, QSplitter, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,QDoubleSpinBox
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
matplotlib.use('Qt5Agg')
import numpy as np
from numpy import linspace

from flim_labs_api import FlimLabsApi

class MplCanvas(FigureCanvas):
    def __init__(self, parent=None, width=5, height=4, dpi=100, title='', nrows=1, ncols=1):
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = []
        for i in range(nrows * ncols):
            self.axes.append(fig.add_subplot(nrows, ncols, i + 1))
        fig.suptitle(title)
        super(MplCanvas, self).__init__(fig)
        
        


class MainWindow(QMainWindow):
    def __init__(self, *args, **kwargs):
        super(MainWindow, self).__init__(*args, **kwargs)

        self.laser_mhz = 40
        self.laser_period_in_nanoseconds = 1000 / self.laser_mhz
        self.acquisition_time_in_seconds_ref = 20
        self.acquisition_time_in_seconds_data = 20
        self.refresh_time_seconds = 5
        self.harmonic_value =  1
        self.tau_phase_reference = 4
        self.tau_modulation_reference = 4 
        self.k = 1 / (2*np.pi*self.harmonic_value*self.laser_mhz*1e6)
        self.m = np.sqrt(1 / (1 + ((self.tau_modulation_reference*1e-9/self.k)**2)))    
        self.phi = math.atan2(self.tau_modulation_reference*1e-9,self.k)
        
                
        self.y_data_ref = np.ones(256)        
        self.x_data = linspace(0, 1000 / self.laser_mhz, 256)
        self.y_data = np.ones(256)
        self.y_data_upd = np.ones(256)
        
        
        ## TO BE VERIFIED
        self.sine_reference = np.sin(2*np.pi*self.harmonic_value*self.x_data / self.laser_period_in_nanoseconds) ## AR
        self.cosine_reference = np.cos(2*np.pi*self.harmonic_value*self.x_data / self.laser_period_in_nanoseconds) ## AR
        #self.sine_reference = np.sin(2*np.pi*self.harmonic_value*self.x_data / 256) 
        #self.cosine_reference = np.cos(2*np.pi*self.harmonic_value*self.x_data / 256) 
        self.s_reference = (self.sine_reference*self.y_data_ref).sum()/self.y_data_ref.sum() 
        self.g_reference = (self.cosine_reference*self.y_data_ref).sum()/self.y_data_ref.sum() 
        
        self.phi_reference = math.atan2(self.s_reference,self.g_reference)
        self.m_reference = np.sqrt(self.g_reference**2+self.s_reference**2)
        self.phi_instr = self.phi_reference - self.phi
        self.m_instr = self.m_reference/self.m 
        
        self.setGeometry(512, 512, 1920, 1080)
        self.setWindowTitle('Phasor analysis ')
        
        #self.laser_measure = False
        
        #if self.laser_measure:
         #  self.api = FlimLabsApi()
          # self.api.set_consumer_handler(self.receive_measure)
        #else:
        self.api = FlimLabsApi()
        self.api.set_consumer_handler(self.receive_point)
        
        
        self.acquiring_ref = False
        self.acquiring_data = False
        
        
        self.tau_life = 3.134 
        
        #list for g and s points to be plotted all together in canvas3.ax2
        self.g_data_referenced_list = []
        self.s_data_referenced_list = []
        
        # Initialize lists to store acquired data
        #self.x_data_list = []
        #self.y_data_list = []
        self.g_data_excel_list = []
        self.s_data_excel_list = []
        self.m_fluo_list = []
        self.phi_fluo_list = []
        
        #initialize the counter for the batches of data in canvas2.ax2
        self.batch_counter = 0

        


        # Create three MplCanvas objects
        self.canvas1 = MplCanvas(self, width=5, height=4, dpi=100, nrows=1, ncols=2)  
        self.canvas1.ax1 = self.canvas1.axes[0]
        self.canvas1.ax1.set_visible(False)        
        self.canvas1.ax2 = self.canvas1.axes[1]
        self.canvas1.ax2.plot(self.x_data, self.y_data_ref) 
        self.canvas1.ax2.set_title('TCSPC calibration sample')        
        self.canvas1.ax2.set_ylabel('Counts')
        self.canvas1.setFixedSize(1920, 330)
                
        self.canvas2 = MplCanvas(self, width=5, height=4, dpi=100, nrows=1, ncols=2)
        self.canvas2.ax1 = self.canvas2.axes[0]
        self.canvas2.ax1.plot(self.x_data, self.y_data)
        self.canvas2.ax1.set_title('TCSPC data sample')
        self.canvas2.ax2 = self.canvas2.axes[1]
        self.canvas2.ax2.plot(self.x_data, self.y_data_upd)
        self.canvas2.ax2.set_title('TCSPC data sample')
        self.canvas2.setFixedSize(1920, 330)

        self.canvas3 = MplCanvas(self, width=5, height=4, dpi=100, nrows = 1, ncols = 2)
        self.canvas3.ax1 = self.canvas3.axes[0]
        #self.canvas3.ax1 = self.canvas3.figure.add_axes([0.1, 0.3, 0.35, self.subplot_height])
        self.canvas3.ax1.plot(self.x_data, self.y_data)
        self.canvas3.ax1.set_xlabel('G')
        self.canvas3.ax1.set_ylabel('S')
        self.canvas3.ax1.set_xlim([-0.005,1.2])
        self.canvas3.ax1.set_ylim([0, 0.6])
        self.canvas3.ax1.set_title('Global phasor plot')
        self.x_circle = np.linspace(-1.0, 1.0, 100)
        self.y_circle = np.linspace(-1.0, 1.0, 100)
        self.X, self.Y = np.meshgrid(self.x_circle,self.y_circle)
        self.F = (self.X**2 + self.Y**2 - self.X)
        self.canvas3.ax1.contour(self.X,self.Y,self.F,[0],colors='b',linewidths=3)
        self.canvas3.ax2 = self.canvas3.axes[1]
        #self.canvas3.ax2 = self.canvas3.figure.add_axes([0.55, 0.3, 0.35, self.subplot_height])
        self.canvas3.ax2.plot(self.x_data, self.y_data)
        self.canvas3.ax2.set_xlabel('G')
        self.canvas3.ax2.set_ylabel('S')
        self.canvas3.ax2.set_xlim([-0.005,1.2])
        self.canvas3.ax2.set_ylim([0, 0.6])
        self.canvas3.ax2.set_title('Global phasor plot')
        self.x_circle_2 = np.linspace(-1.0, 1.0, 100)
        self.y_circle_2 = np.linspace(-1.0, 1.0, 100)
        self.X_2, self.Y_2 = np.meshgrid(self.x_circle_2,self.y_circle_2)
        self.F_2 = (self.X_2**2 + self.Y_2**2 - self.X_2)        
        self.canvas3.ax2.contour(self.X_2,self.Y_2,self.F_2,[0],colors='b',linewidths=3)
        self.canvas3.setFixedSize(1920, 330)
        
        # Create a vertical layout and add the three canvas widgets to it
        layout = QVBoxLayout()
        layout.addWidget(self.canvas1)
        layout.addWidget(self.canvas2)
        layout.addWidget(self.canvas3)

        # Create a central widget and set the layout
        central_widget = QWidget()
        central_widget.setLayout(layout)
        self.setCentralWidget(central_widget)

        # draw a button to start the acquisition in reference
        self.start_button_ref = QPushButton('Start', self)
        self.start_button_ref.move(1050, 5)
        self.start_button_ref.clicked.connect(self.start_acquisition)
        
        # draw a button to STOP the acquisition in reference
        self.stop_button_ref = QPushButton('Stop', self)
        self.stop_button_ref.move(1175, 5)
        self.stop_button_ref.clicked.connect(self.stop_acquisition)
        
        # draw a button to start the acquisition in data
        self.start_button_data = QPushButton('Start', self)
        self.start_button_data.move(5, 330)
        self.start_button_data.clicked.connect(self.start_acquisition)
        
        
        # draw a button to stop the acquisition in data
        self.stop_button_data = QPushButton('Stop', self)
        self.stop_button_data.move(120, 330)
        self.stop_button_data.clicked.connect(self.stop_acquisition)
        self.stop_button_data.clicked.connect(self.save_to_excel)
        self.mutex = QMutex()
        
        
        # create a label and spin box for laser frequency setting (sync out)
        self.freq_label = QLabel('Laser frequency(MHz):', self)
        self.freq_label.move(200, 45)
        self.freq_label.adjustSize()

        self.freq_spinbox = QSpinBox(self)
        self.freq_spinbox.move(200, 70)
        self.freq_spinbox.setMinimum(20)
        self.freq_spinbox.setMaximum(80)
        self.freq_spinbox.setValue(self.laser_mhz)
        self.freq_spinbox.valueChanged.connect(self.set_laser_frequency)
        
        #create a button to measure the laser frequency with measure-frequency method
        self.sync_laser_in = QPushButton('Sync Laser', self)
        self.sync_laser_in.move(400,45)
        self.sync_laser_in.clicked.connect(self.frequency_meter)
        
        #create a label to display the result of the laser frequency measurement 
        self.measure_label = QLabel(self)
        self.measure_label.move(310, 145)
        self.measure_label.setText('Press Sync Laser to measure Laser frequency')
        self.measure_label.adjustSize()

        
        # create a label and spin box for acquisition time setting in reference
        self.time_label = QLabel('Acquisition Time (Seconds):', self)
        self.time_label.move(5, 45)
        self.time_label.adjustSize()

        self.time_spinbox = QSpinBox(self)
        self.time_spinbox.move(5, 70)
        self.time_spinbox.setMinimum(1)
        self.time_spinbox.setMaximum(300)
        self.time_spinbox.setValue(self.acquisition_time_in_seconds_ref)
        self.time_spinbox.valueChanged.connect(self.set_acquisition_time_ref)
        
        # create a label and spin box for acquisition time setting in data
        self.time_label = QLabel('Acquisition Time (Seconds):', self)
        self.time_label.move(5, 370)
        self.time_label.adjustSize()

        self.time_spinbox = QSpinBox(self)
        self.time_spinbox.move(5, 395)
        self.time_spinbox.setMinimum(1)
        self.time_spinbox.setMaximum(300)
        self.time_spinbox.setValue(self.acquisition_time_in_seconds_data)
        self.time_spinbox.valueChanged.connect(self.set_acquisition_time_data)
        
        # create a label to show the photons received in ref
        self.phase_label_ref = QLabel(self)
        self.phase_label_ref.move(5, 245 )
        self.phase_label_ref.setText('Points received: 0')
        self.phase_label_ref.adjustSize()
        self.points_received_ref = 0
        
        # create a label to show the photons received in data
        self.phase_label = QLabel(self)
        self.phase_label.move(5, 550)
        self.phase_label.setText('Points received: 0')
        self.phase_label.adjustSize()
        self.points_received = 0
        
        #create a label and spin box to set the refresh time of the histogram in data
        self.refresh_label = QLabel('Refresh Time (Seconds):', self)
        self.refresh_label.move(5, 480)
        self.refresh_label.adjustSize()
        
        self.refresh_spinbox = QSpinBox(self)
        self.refresh_spinbox.move(5, 505)
        self.refresh_spinbox.setMinimum(3)
        self.refresh_spinbox.setMaximum(30)
        self.refresh_spinbox.setValue(self.refresh_time_seconds )
        self.refresh_spinbox.valueChanged.connect(self.set_refresh_time_in_seconds)
        
        # create timer to refresh histogram of reference
        self.refresh_timer = QTimer()
        self.refresh_timer.timeout.connect(self.refresh_histogram_ref)
        self.refresh_timer.setInterval(100)
        self.refresh_timer.start()
        
        # create timer to refresh histogram of data
        self.refresh_timer_data = QTimer()
        self.refresh_timer_data.timeout.connect(self.refresh_histogram)
        self.refresh_timer_data.setInterval(100)
        self.refresh_timer_data.start()
        
        # create timer to refresh histogram of data refreshed
        #self.refresh_timer_data_refreshed = QTimer()
        #self.refresh_timer_data_refreshed.timeout.connect(self.refresh_histogram_refreshed)
        #self.refresh_timer_data_refreshed.setInterval(100)
        #self.refresh_timer_data_refreshed.start()
        
        # create a timer to update canvas2 with new data every tot seconds 
        self.timer = QTimer(self)
        self.timer.timeout.connect(self.update_canvas2)
        self.timer.setInterval(self.refresh_time_seconds * 1000)
        self.timer.start()
        
        #draw a panel for harmonic selection in phasor analysis
        self.harmonic_label = QLabel('Harmonic:', self)
        self.harmonic_label.move(5,125)
        self.harmonic_label.adjustSize()

        self.harmonic_spinbox = QSpinBox(self)
        self.harmonic_spinbox.move(5, 145)
        self.harmonic_spinbox.setMinimum(1)
        self.harmonic_spinbox.setMaximum(10)
        self.harmonic_spinbox.setValue(self.harmonic_value)
        self.harmonic_spinbox.valueChanged.connect(self.set_harmonic_value)
        
        #draw a panel to show the tau_phase_reference
        self.tau_label = QLabel('tau_phase (ns):', self)
        self.tau_label.move(600, 45)
        self.tau_label.adjustSize()

        self.tau_spinbox = QDoubleSpinBox(self)
        self.tau_spinbox.move(600, 70)
        self.tau_spinbox.setMinimum(1.0)
        self.tau_spinbox.setMaximum(10.0)
        self.tau_spinbox.setSingleStep(0.01)  # Set the step size for decimal values
        self.tau_spinbox.setValue(self.tau_phase_reference)
        self.tau_spinbox.valueChanged.connect(self.set_tau_phase_reference)

        
        #draw a panel to show the tau_modulation_reference
        self.tau_mod_label = QLabel('tau_modulus (ns):', self)
        self.tau_mod_label.move(800, 45)
        self.tau_mod_label.adjustSize()

        self.tau_mod_spinbox = QDoubleSpinBox(self)
        self.tau_mod_spinbox.move(800, 70)
        self.tau_mod_spinbox.setMinimum(1)
        self.tau_mod_spinbox.setMaximum(10)
        self.tau_mod_spinbox.setSingleStep(0.01) 
        self.tau_mod_spinbox.setValue(self.tau_modulation_reference )
        self.tau_mod_spinbox.valueChanged.connect(self.set_tau_modulus_reference)
        
        # create a label to show the phase of expected sample
        self.tau_phase_label_ref = QLabel(self)
        self.tau_phase_label_ref.move(600, 145)
        self.tau_phase_label_ref.setText('phase =' +str(np.round(self.phi,4)))
        self.tau_phase_label_ref.adjustSize()
        
        # create a label to show the modulus of expected sample
        self.tau_modulus_label_ref = QLabel(self)
        self.tau_modulus_label_ref.move(800, 145)
        self.tau_modulus_label_ref.setText('modulus =' +str(np.round(self.m,4)))
        self.tau_modulus_label_ref.adjustSize()
        
        # create a label to show the modulus of the reference
        self.tau_ref_modulus_label_ref = QLabel(self)
        self.tau_ref_modulus_label_ref.move(200, 245)
        self.tau_ref_modulus_label_ref.setText('IRF modulus =' +str(np.round(self.m_instr,4)))
        self.tau_ref_modulus_label_ref.adjustSize()
        
        # create a label to show the phase of the reference
        self.tau_ref_phase_label_ref = QLabel(self)
        self.tau_ref_phase_label_ref.move(400, 245)
        self.tau_ref_phase_label_ref.setText('IRF phase =' +str(np.round(self.phi_instr,4)))
        self.tau_ref_phase_label_ref.adjustSize()
        #create a timer to refresh the phasor plot
        #self.refresh_timer_phasor = QTimer()
        #self.refresh_timer_phasor.timeout.connect(self.refresh_phasor_plot)
        #self.refresh_timer_phasor.setInterval(6000)
        #self.refresh_timer_phasor.start()
        
        # create a label to show the g value
        self.g_label = QLabel(self)
        self.g_label.move(5,760)
        self.g_label.setText('g = 0')
        self.g_label.adjustSize()
        
        # create a label to show the s value
        self.s_label = QLabel(self)
        self.s_label.move(5, 810)
        self.s_label.setText('s = 0')
        self.s_label.adjustSize()
        
        # create a label to show the tau value of lifetime
        self.tau_life_phase_label = QLabel(self)
        self.tau_life_phase_label.move(5, 860)
        self.tau_life_phase_label.setText('tau_phase = 0')
        #self.tau_life_phase_label.setText('tau_phase =' +str(np.round(self.tau_life*1e9,3)))
        self.tau_life_phase_label.setText('tau_phase = 0')
        self.tau_life_phase_label.adjustSize()
        
        # create a label to show the tau value of lifetime
        self.tau_life_mod_label = QLabel(self)
        self.tau_life_mod_label.move(5, 910)
        self.tau_life_mod_label.setText('tau_modulation =' +str(np.round(self.tau_life*1e9,3)))
        self.tau_life_mod_label.adjustSize()
        
        
        
        self.show()
        
    def closeEvent(self, event):
        np.savetxt('Reference_TCSPC.txt', np.vstack((self.x_data, self.y_data_ref)).T)
        np.savetxt('Data_TCSPC.txt', np.vstack((self.x_data, self.y_data)).T)
        self.api.stop_acquisition()
        event.accept()
        
    def set_laser_frequency(self, value):
        self.laser_mhz = value  
        self.x_data = linspace(0, 1000 / self.laser_mhz, 256)
        self.laser_period_in_nanoseconds = 1000 / self.laser_mhz
        
    

                        
    def set_acquisition_time_ref(self, value):
        self.acquisition_time_in_seconds_ref = value

    def set_acquisition_time_data(self, value):
        self.acquisition_time_in_seconds_data = value  
        self.w = linspace(0, self.acquisition_time_in_seconds_data / self.acquisition_time_in_seconds_data, 256) 
        
    def set_refresh_time_in_seconds(self, value):
        self.refresh_time_seconds = value
        self.canvas2.ax2.set_title(f'{self.refresh_time_seconds} seconds TCSPC data sample')
        self.canvas3.ax2.set_title(f'{self.refresh_time_seconds} seconds Refreshed phasor plot')

        
    
    def set_harmonic_value(self, value):
        self.harmonic_value = value
        self.k = 1/(2*np.pi*self.harmonic_value*self.laser_mhz*1e6)
        self.phi = math.atan2(self.tau_phase_reference*1e-9,self.k)
        self.tau_phase_label_ref.setText(f'phase : {self.phi:,}')
        self.m = np.sqrt(1/(1+((self.tau_modulation_reference*1e-9/self.k)**2)))        
        self.tau_modulus_label_ref.setText(f'modulus : {self.m:,}')
        
    def set_tau_phase_reference(self, value):
        self.tau_phase_reference = value 
        self.phi = math.atan2(self.tau_phase_reference*1e-9,self.k)
        self.tau_phase_label_ref.setText(f'phase : {self.phi:,}')        
        
    def set_tau_modulus_reference(self, value):
        self.tau_modulation_reference = value   
        self.m = np.sqrt(1/(1+((self.tau_modulation_reference*1e-9/self.k)**2)))        
        self.tau_modulus_label_ref.setText(f'modulus : {self.m:,}')
        
    def start_acquisition(self):
        sender = self.sender()
        if sender == self.start_button_ref:
            self.start_button_ref.setEnabled(False)
            self.acquiring_ref = True
            self.y_data_ref = np.zeros(256)
            self.api.set_firmware("firmwares\\spectroscopy_40MHz.flim")
            self.api.acquire_spectroscopy(
            laser_frequency_mhz=self.laser_mhz,
            acquisition_time_seconds=self.acquisition_time_in_seconds_ref
            )
            self.stop_button_ref.setEnabled(True)
        elif sender == self.start_button_data:
            self.start_button_data.setEnabled(False)
            self.acquiring_data = True
            self.y_data = np.zeros(256)
            self.api.set_firmware("firmwares\\spectroscopy_40MHz.flim")
            self.api.acquire_spectroscopy(
            laser_frequency_mhz=self.laser_mhz,
            acquisition_time_seconds=self.acquisition_time_in_seconds_data
            )
            self.stop_button_data.setEnabled(True)

    
    def stop_acquisition(self):
        sender = self.sender()
        if sender == self.stop_button_ref:
            self.stop_button_ref.setEnabled(False)
            self.acquiring_ref = False
            self.api.stop_acquisition()
            #self.y_data_ref = np.zeros(256)
            #self.points_received_ref = 0
            self.start_button_ref.setEnabled(True)
        elif sender == self.stop_button_data:    
            self.stop_button_data.setEnabled(False)
            self.acquiring_data = False
            self.api.stop_acquisition()
            #self.y_data = np.zeros(256)
            self.points_received = 0
            self.start_button_data.setEnabled(True) 
        
        
    def receive_point(self, channel, time_bin, micro_time, monotonic_counter, macro_time):
        if self.acquiring_ref:
            self.mutex.lock()
            self.y_data_ref[(time_bin+90)%256] += 1
            self.points_received_ref += 1            
            self.mutex.unlock()
        elif self.acquiring_data:
            self.mutex.lock()
            self.y_data[(time_bin+90)%256] += 1
            self.y_data_upd[(time_bin+90)%256] += 1
            self.points_received += 1
            self.mutex.unlock()   
            
            #self.x_data_list.append(self.x_data)
            #self.y_data_list.append(self.y_data)    
            
    def frequency_meter(self):
        self.api.set_consumer_handler(self.receive_measure)
        self.measure_label.setText('Waiting for measure...')
        self.measure_label.adjustSize()
        self.sync_laser_in.setEnabled(False)     
        self.api.set_firmware("firmwares\\frequency_meter.flim")
        self.api.acquire_measure_frequency()
        self.update()        
        
    def receive_measure(self, frequency):
        self.measure_label.setText('Measurement: ' + str(frequency) + ' MHz')
        self.measure_label.adjustSize()
        self.sync_laser_in.setEnabled(True)
        self.api.stop_acquisition()
        self.api.set_consumer_handler(self.receive_point)
        self.update()    
    
               

    def refresh_histogram_ref(self):
        self.canvas1.ax2.clear()
        self.canvas1.ax2.plot(self.x_data, self.y_data_ref)
        #self.chart.axes.set_xlabel('Time (ns)')
        #self.chart.axes.set_ylabel('Photon counts')
        self.canvas1.ax2.set_xlim([0, 1000 / self.laser_mhz])
        self.canvas1.ax2.set_title('TCSPC calibration sample')
        self.canvas1.ax2.set_xlabel('Time (ns)')
        self.canvas1.ax2.set_ylabel('Counts')
        self.canvas1.draw()
        # format points received with commas
        self.phase_label_ref.setText(f'Total photons: {self.points_received_ref:,}')
        self.s_reference = (self.sine_reference*self.y_data_ref).sum()/self.y_data_ref.sum() 
        self.g_reference = (self.cosine_reference*self.y_data_ref).sum()/self.y_data_ref.sum() 
        
        self.phi_reference = math.atan2(self.s_reference,self.g_reference)
        self.m_reference = np.sqrt(self.g_reference**2+self.s_reference**2)
        self.phi_instr = self.phi_reference - self.phi
        self.m_instr = self.m_reference/self.m 
        self.tau_ref_modulus_label_ref.setText('IRF modulus =' +str(np.round(self.m_instr,4)))
        self.tau_ref_phase_label_ref.setText('IRF phase =' +str(np.round(self.phi_instr,4)))
        self.phase_label_ref.adjustSize()
        
    def refresh_histogram(self):
        self.canvas2.ax1.clear()
        self.canvas2.ax1.plot(self.x_data, self.y_data)
        #self.chart.axes.set_xlabel('Time (ns)')
        #self.chart.axes.set_ylabel('Photon counts')
        self.canvas2.ax1.set_xlim([0, 1000 / self.laser_mhz])
        self.canvas2.ax1.set_title('Global decay histogram')
        self.canvas2.ax1.set_xlabel('Time (ns)')
        self.canvas2.ax1.set_ylabel('Counts')
        self.canvas2.draw_idle()
        # format points received with commas
        self.phase_label.setText(f'Total photons: {self.points_received:,}')
        self.phase_label.adjustSize() 
        
     
    def update_canvas2(self):
        self.mutex.lock()
        self.canvas2.ax2.clear()
        self.canvas2.ax2.plot(self.x_data, self.y_data_upd)
        self.canvas2.ax2.set_xlim([0, 1000 / self.laser_mhz])
        self.canvas2.ax2.set_title(f'{self.refresh_time_seconds} seconds TCSPC data sample')
        self.canvas2.ax2.set_xlabel('Time (ns)')
        self.canvas2.ax2.set_ylabel('Counts')
        self.canvas2.draw_idle()
        
        if self.acquiring_data:
           QTimer.singleShot(self.refresh_time_seconds * 1000, self.save_and_reset_data)
           
        S = np.sum(self.y_data_upd * self.sine_reference) / np.sum(self.y_data_upd)
        G = np.sum(self.y_data_upd * self.cosine_reference) / np.sum(self.y_data_upd)
        m_data = np.sqrt(G**2+S**2) 
        phi_data = np.arctan2(S,G) 
         
        m_fluo = m_data/self.m_instr  
        phi_fluo = phi_data - self.phi_instr
        g_data_referenced = m_fluo*np.cos(phi_fluo) 
        s_data_referenced = m_fluo*np.sin(phi_fluo)  
        
        if self.acquiring_data:
           self.g_data_excel_list.append(g_data_referenced)
           self.s_data_excel_list.append(s_data_referenced)
           self.m_fluo_list.append(m_fluo)
           self.phi_fluo_list.append(phi_fluo)
              
        #print(len(self.g_data_excel_list)) 
        #print(len(self.s_data_excel_list))     
        #print(len(self.m_fluo_list)) 
        #print(len(self.m_fluo_list))         
        
        S_1 = np.sum(self.y_data * self.sine_reference) / np.sum(self.y_data)
        G_1 = np.sum(self.y_data * self.cosine_reference) / np.sum(self.y_data)
        m_data_1 = np.sqrt(G_1**2 + S_1**2) 
        phi_data_1 = np.arctan2(S_1,G_1) 
         
        m_fluo_1 = m_data_1 / self.m_instr  
        phi_fluo_1 = phi_data_1 - self.phi_instr
        g_data_referenced_1 = m_fluo_1 * np.cos(phi_fluo_1) 
        s_data_referenced_1 = m_fluo_1 * np.sin(phi_fluo_1) 
        self.g_label.setText('g =' +str(np.round(g_data_referenced_1,3)))
        self.s_label.setText('s =' +str(np.round(s_data_referenced_1,3)))
        self.g_label.adjustSize()
        self.s_label.adjustSize()
        
        
        phi_fluo_avg = np.arctan2(s_data_referenced_1, g_data_referenced_1)
        m_fluo_avg = np.sqrt(g_data_referenced_1**2 + s_data_referenced_1**2)

        tau_p = self.k*math.tan(phi_fluo_avg)
        
        par = (1/m_fluo_avg**2)
        tau_m = self.k*np.sqrt(par-1)
        
        self.tau_life_phase_label.setText('tau_phase =' + str(np.round(tau_p*1e9,3)))
        self.tau_life_mod_label.setText('tau_modulation =' + str(np.round(tau_m*1e9,3)))
        
        if self.acquiring_data:
           self.g_data_referenced_list.append(g_data_referenced)
           self.s_data_referenced_list.append(s_data_referenced)
        
        # Update the phasor plot di sinistra (un solo punto)
        self.canvas3.ax1.clear()
        self.canvas3.ax1.set_xlim([-0.005,1.2])
        self.canvas3.ax1.set_ylim([0, 0.6])
        self.canvas3.ax1.plot(g_data_referenced_1, s_data_referenced_1, 'bo')  #non mi mostrava il plot perchè a denominatore G e S hanno 0  non avendo ancora acquisito i dati e avendoli inizializzati a zero
        #self.canvas3.ax1.set_xlim([0, 1])
        self.canvas3.ax1.set_title('Global phasor plot')
        self.canvas3.ax1.set_xlabel('g')
        self.canvas3.ax1.set_ylabel('s')
        self.canvas3.ax1.contour(self.X_2,self.Y_2,self.F_2,[0],colors='b',linewidths=3)
        # Update the phasor plot di destra
        
        self.canvas3.ax2.clear()
        self.canvas3.ax2.set_xlim([-0.005, 1.2])
        self.canvas3.ax2.set_ylim([0, 0.6])
        self.canvas3.ax2.plot(self.g_data_referenced_list, self.s_data_referenced_list, 'bo')
        self.canvas3.ax2.set_title(f'{self.refresh_time_seconds} seconds Refreshed phasor plot')
        self.canvas3.ax2.set_xlabel('g')
        self.canvas3.ax2.set_ylabel('s')
        self.canvas3.ax2.contour(self.X, self.Y, self.F, [0], colors='b', linewidths=3)
        
        
       
        self.canvas3.draw()
        
        
        self.mutex.unlock()

        
    
    def save_and_reset_data(self):
        filename = f'data_batch_{self.batch_counter}.txt'
        np.savetxt(filename, np.vstack((self.x_data, self.y_data_upd)).T)
        self.batch_counter += 1
        self.y_data_upd = np.zeros_like(self.y_data_upd)  #qui è la chiave del problema   
        self.canvas2.ax2.clear()
        self.canvas2.ax2.plot(self.x_data, self.y_data_upd)
        self.canvas2.ax2.set_title('TCSPC')
        self.canvas2.ax2.set_xlabel('Time (ns)')
        self.canvas2.ax2.set_ylabel('Counts')
        self.canvas2.draw_idle()
        #QTimer.singleShot(self.refresh_time_seconds * 1000, self.update_canvas2) #prova a moltiplicare refreshed_time per 1000
        
    def save_to_excel(self):
        #for _ in range(len(self.x_data_list) - len(self.g_data_excel_list)):
           
         #   self.g_data_excel_list.append(0)
          #  self.s_data_excel_list.append(0)
           # self.m_fluo_list.append(0)
            #self.phi_fluo_list.append(0)
        data = {
            #'time bin': self.x_data_list,
            #'photon count': self.y_data_list,
            'g_data': self.g_data_excel_list,
            's_data': self.s_data_excel_list,
            'm_fluo': self.m_fluo_list,
            'phi_fluo': self.phi_fluo_list
        }
        
        #data_1 = {
        
         #   'time bin': self.x_data_list,
          #  'photon count': self.y_data_list
        #}

        df = pd.DataFrame(data)
        #df_1 = pd.DataFrame(data_1)
        df.to_excel('phasors_data.xlsx', index=False)    
        #df_1.to_excel('TCSPC_data.xlsx', index=False)
        
    
                




        
        

if __name__ == '__main__':
    app = QApplication(sys.argv)
    window = MainWindow()
    window.show()
    sys.exit(app.exec_())