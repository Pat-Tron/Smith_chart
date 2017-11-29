# -*- coding:utf-8 -*-

import numpy as np
import matplotlib.pyplot as plt

PI = 3.14159

class Smith_co:
    """
    Simple Smith Chart coordinate drawing.

    User guide:

    After opening application ,input characteristic impedance in integer and
    impedance_to_be_matched in complex both with unit ohm, then input frequency
    with unit Mhz.Press enter directly will use deflaut value :

        deflaut Characteristic Impedance : 50 Ohms
        deflaut Impedance_to_be_matched : 50 + 0j Ohms
        deflaut Frequency : 1 MHz

    Chr_impedance and frequency can't change after you press the 3rd enter.Only
    impedance_to_be_matched can be redefine with mouse clicking and it is
    unaccuratly but conveniently.

    Values and chart label will display on the right. Mouse's (x, y) postion can
    be found on the leftdown corner (also means reflection coefficient).

    (Switch mouse to pan could activate click motion.........)

    """

    def __init__(self):
        """
        Creates the figure and draws the grid.
        """


        # Welcome and intial data input.
        print("\n Thanks to use SmithChart! \n")

        z0_t = input("Please input characteristic impedance (integer):")
        self.z0 = complex(z0_t if z0_t else 50)

        im_t = input("And impedance to be matched (complex):")
        self.impedance2match = complex(im_t if im_t else self.z0)

        fr_t = input("And frequency of system (Unit : MHz):")
        self.frequency = 10 ** 6 * float(fr_t if fr_t else 1)
        self.omega = self.frequency * 2 * PI


        # Basic chart initialization.
        self.fig = plt.figure(figsize = (17, 12.5))
        self.fig.subplots_adjust(left=0.1, bottom=0.1, right=0.7, top=0.9)
        plt.title('Smith Chart', fontsize = 25)
        plt.xlabel('Real part of reflection coefficient', fontsize = 22)
        plt.ylabel('Imaginary part of reflection coefficient', fontsize = 22)
        plt.xticks(fontsize = 16)
        plt.yticks(fontsize = 16)

        # Chart label and value of matching update.
        self.labelsAndText()

        # Initialize impedance matching path lines.
        self.lineInit()

        # Smith round coordinate 's drawing.
        self.drawGrid()


        # Path drawing of according to input value.
        # Draw input impedance point, green.
        self.z_point = plt.scatter(0, 0, s = 150, c ='#165900')

        # 1st order impedance matching.
        if self.impedance2match.real == 0:
            y = self.impedance2match.imag / self.z0.real
            self.markZ(self.z2gamma(complex(1, y)))
            self.first_order_matching(self.impedance2match/self.z0)

        # 2nd order impedance matching.
        elif self.impedance2match.imag != 0 and self.impedance2match.real != 1:
            self.markZ(self.z2gamma(self.impedance2match/self.z0))
            self.second_order_matching(self.impedance2match/self.z0)

        # Validate mouse action.
        self.fig.canvas.mpl_connect('button_press_event', self.onclick)

    def labelsAndText(self):
        "Updating of chart label and value of impedance matching."

        c_i = 'Characteristic Impedance : %d Ohms' %self.z0.real
        f = 'Frequency : %d Mhz' % (self.frequency / (10 ** 6))
        abbr = 's - Series     p - Parallel'

        self.i_text = plt.text(1.15, 1,   'Impedance  : 0 + 0 j', fontsize = 19)
        self.a_text = plt.text(1.15, 0.92,'Admittance : 0 + 0 j', fontsize = 19)
        self.g_text = plt.text(1.15, 0.84,'Gamma      : 0 + 0 j', fontsize = 19)
        self.c_text = plt.text(1.15, 0.65, c_i, fontsize = 16)
        self.f_text = plt.text(1.15, 0.59, f, fontsize = 16)
        self.b_text = plt.text(1.15, 0.53, abbr, fontsize = 16, weight = 'bold')

        self.text_1 = plt.text(1.15, 0.3, '1st order matching :', fontsize = 19)
        self.f_rea_t = plt.text(1.15, 0,'____ Reactance : 0j Ohms',
            color = '#E5C000', fontsize = 17, weight = 'bold')

        self.text_2 = plt.text(1.15, -0.2, '2nd order matching :',fontsize = 19)
        self.s_l_2_t = plt.text(1.15, -0.3,'-.-.- switch to admittance chart',
            color = '#E06900', fontsize = 15)
        self.s_l_3_t = plt.text(1.15, -0.37,'-.-.- switch to another chart',
            color = '#00BF9C', fontsize = 15)
        self.s_l_4_t = plt.text(1.15, -0.44,'-.-.- switch to another chart',
            color = '#0061C2', fontsize = 15)
        self.s_r_1_t = plt.text(1.15, -0.64,'----- Reactance_1_1 : 0j Ohms',
            color = '#E08B00', fontsize = 17, weight = 'bold')
        self.s_r_2_t = plt.text(1.15, -0.74,'----- Reactance_1_2 : 0j Ohms',
            color = '#0061C2', fontsize = 17, weight = 'bold')
        self.s_r_3_t = plt.text(1.15, -0.94,'----- Reactance_2_1 : 0j Ohms',
            color = '#F03F00', fontsize = 17, weight = 'bold')
        self.s_r_4_t = plt.text(1.15, -1.04,'----- Reactance_2_2 : 0j Ohms',
            color = '#00BF9C', fontsize = 17, weight = 'bold')

        # Every axis/circle's label.
        self.circleLable('0', 0, 0, -0.05, 0.01)
        self.circleLable('0.1', 0.1, 0, -0.1, 0.01)
        self.circleLable('0.3', 0.3, 0, -0.1, 0.01)
        self.circleLable('0.6', 0.6, 0, -0.1, 0.01)
        self.circleLable('1', 1, 0, -0.05, 0.01)
        self.circleLable('1.5', 1.5, 0, -0.1, 0.01)
        self.circleLable('2.4', 2.4, 0, -0.1, 0.01)
        self.circleLable('4.4', 4.4, 0, -0.1, 0.01)
        self.circleLable('10', 10, 0, 0.1, 0.01)
        self.circleLable('0.3j', 0, 0.3, -0.15, 0)
        self.circleLable('0.6j', 0, 0.6, -0.15, 0)
        self.circleLable('j', 0, 1, 0, 0.02)
        self.circleLable('1.5j', 0, 1.5, 0, 0.01)
        self.circleLable('2.5j', 0, 2.5, 0, 0.01)
        self.circleLable('6j', 0, 6, 0, 0)
        self.circleLable('-0.3j', 0, -0.3, -0.15, -0.05)
        self.circleLable('-0.6j', 0, -0.6, -0.15, -0.05)
        self.circleLable('-j', 0, -1, -0.02, -0.06)
        self.circleLable('-1.5j', 0, -1.5, 0, -0.05)
        self.circleLable('-2.5j', 0, -2.5, 0, -0.05)
        self.circleLable('-6j', 0, -6, 0, -0.02)
        self.circleLable('inf', 0, 0, 2, 0)

    def lineInit(self):
        "Initialize impedance matching path lines."

        # 1st order matching, reactance, yellow
        self.f_order_reactance, = plt.plot([], [], c = '#FFDE02', lw = 7)

        # 2nd order matching, no.1 reactance with positive imaginary part,orange
        self.s_order_reactance_1_1,=plt.plot([],[],c='#E08B00', lw = 5, ls='--')
        # 2nd order matching, no.1 reactance with negative imaginary part, red
        self.s_order_reactance_1_2,=plt.plot([],[],c='#F03F00', lw = 5, ls='--')
        # 2nd order matching, no.2 reactance with positive imaginary part, cyan
        self.s_order_reactance_2_1,=plt.plot([],[],c='#00BF9C', lw = 5, ls='--')
        # 2nd order matching, no.2 reactance with negative imaginary part, blue
        self.s_order_reactance_2_2,=plt.plot([],[],c='#0061C2', lw = 5, ls='--')
        # 2nd order matching, swith line connected to no.1 reactance with
        # negative imaginary part, blue
        self.s_order_switch_1,= plt.plot([], [], c = '#0061C2', lw = 2, ls='-.')
        # 2nd order matching, swith line connected to no.1 reactance with
        # positive imaginary part, cyan
        self.s_order_switch_2,= plt.plot([], [], c = '#00BF9C', lw = 2, ls='-.')
        # 2nd order matching, swith line may or maynot appear, red
        self.s_order_switch_t,= plt.plot([], [], c = '#E06900', lw = 2, ls='-.')

    def drawGrid(self):
        """
        Draws the Smith Chart grid.
        """
        resistance = [
        0.02, 0.04, 0.06, 0.08,
        0.12, 0.15, 0.18, 0.21, 0.24, 0.27,
        0.35, 0.40, 0.45, 0.50, 0.55,
        0.68, 0.76, 0.84, 0.92,
        1.1, 1.2, 1.3, 1.4,
        1.63, 1.8, 2.0, 2.2,
        2.8, 3.2, 3.6, 4.0,
        5.7, 7, 8.3, 9.6,
        ]

        reactance = [
        0.07, 0.15, 0.23,
        0.37, 0.45, 0.53,
        0.7, 0.8, 0.9,
        1.12, 1.24, 1.36,
        1.75, 2, 2.25,
        3.3, 4, 5,
        10,15
        ]

        self.drawXCircle(0, width = 2)
        self.drawXCircle(0.1)
        self.drawXCircle(0.3)
        self.drawXCircle(0.6)
        self.drawXCircle(1)
        self.drawXCircle(1.5)
        self.drawXCircle(2.4)
        self.drawXCircle(4.4)
        self.drawXCircle(10)
        for i in resistance : self.drawXCircle(i, 'k:')

        self.drawYCircle(0)
        self.drawYCircle(0.3)
        self.drawYCircle(0.6)
        self.drawYCircle(1)
        self.drawYCircle(1.5)
        self.drawYCircle(2.5)
        self.drawYCircle(6)
        self.drawYCircle(-0.3)
        self.drawYCircle(-0.6)
        self.drawYCircle(-1)
        self.drawYCircle(-1.5)
        self.drawYCircle(-2.5)
        self.drawYCircle(-6)
        for i in reactance : self.drawYCircle(i, 'k:')
        for i in reactance : self.drawYCircle(-i, 'k:')

    def onclick(self, event):
        "Response click motion."

        if event.xdata :
            "If in chart square."
            g = complex(event.xdata, event.ydata)
            z = self.gamma2z(g)

            if z.real >= 0:
                "If in chart round."

                if z.real == 0:
                    "Only 1st order matching."
                    self.first_order_matching(self.gamma2z(g))

                elif event.ydata != 0 and z.real != 1:
                    """
                    If impedance is right on real axis or 1 Ohm cicle, then
                    there is no 2nd order matching.
                    """
                    self.second_order_matching(self.gamma2z(g))

                # Update impedance parameter, and mark it.
                self.markZ(g)

    def markZ(self, g):
        "Mark a impedance on the chart."

        z = self.gamma2z(g)
        I = z * self.z0
        A = complex(1/I if I!=0 else 10**6)

        self.z_point.set_offsets((g.real, g.imag))
        self.i_text.set_text('Impedance  : %.4f + %.4f j' %(I.real, I.imag))
        self.a_text.set_text('Admittance : %.4f + %.4f j' %(A.real, A.imag))
        self.g_text.set_text('Gamma      : %.4f + %.4f j' %(g.real, g.imag))

    def first_order_matching(self, z):
        "Draw 1st order impedance(normalized) matching path when it's cascade."

        # Drawing C or L line.
        if z.imag != 0:
            flag = z.imag / abs(z.imag)
            reactance_list = [self.z2gamma(complex(1, flag * y))
                for y in np.logspace(-2, np.log10(abs(z.imag)), 100)]
            self.f_order_reactance.set_data([z.real for z in reactance_list],
                [z.imag for z in reactance_list])

        # Update reactance label.
        if z.imag > 0:
            "Reactance = jwL"
            self.f_rea_t.set_text('____ s_Inductor(L) : %.3f uH' %(
            z.imag / self.omega * (10 ** 6) * self.z0.real))
        elif z.imag < 0:
            "Reactance = 1/jwC"
            self.f_rea_t.set_text('____ s_Capacitor(C) : %.3f nF' %(
            1 / -z.imag / self.omega * (10 ** 9) / self.z0.real))
        else :
            "Reactance = 0"
            self.f_rea_t.set_text('____ Reactance : 0 Ohms')

    def second_order_matching(self, z):
        """
        Draw 2nd order impedance(normalized) matching path.
        Group 1 :orange and blue
        Group 2 :red and cyan
        """

        self.f_order_reactance.set_data([], [])

        if z.real <= 1:
            "Don't need to switch chart, first component is connected in series."

            # First reactance with positive imaginary.
            f_imag = np.sqrt(z.real - z.real ** 2)
            # Second susceptance with positive imaginary.
            s_imag = (-1 / complex(z.real, f_imag)).imag

            #Draw matching path.
            self.second_order_matching_series(z, f_imag)
            self.s_order_switch_t.set_data([],[])

            # Update reactance  and susceptance label.
            # L of group 1 susceptance = 1/jwL
            self.s_r_2_t.set_text('----- p_Inductor(L) : %.3f uH'
            %(1 / s_imag / self.omega * 10 ** 6 * self.z0.real))
            # C of group 2 susceptance = jwC
            self.s_r_4_t.set_text('----- p_Capacitor(C) : %.3f nF'
            %(s_imag / self.omega * 10 ** 9 / self.z0.real))

            if z.imag >= f_imag:
                "Z is upside 1S admittance circle."

                # reactance of group 1 = jwL
                self.s_r_1_t.set_text('____ s_Inductor(L) : %.3f uH'
                %((z.imag - f_imag)/ self.omega * 10 ** 6 * self.z0.real))
                # reactance of group 2 = jwL
                self.s_r_3_t.set_text('____ s_Inductor(L) : %.3f uH'
                %((z.imag + f_imag)/ self.omega * 10 ** 6 * self.z0.real))

            elif z.imag >= -f_imag :
                "Z is inside 1S admittance circle."

                # reactance of group 1 = 1/jwC
                self.s_r_1_t.set_text('____ s_Capacitor(C) : %.3f nF'
                %(1/(f_imag - z.imag)/ self.omega * 10 ** 9 / self.z0.real))
                # reactance of group 2 = jwL
                self.s_r_3_t.set_text('____ s_Inductor(L) : %.3f uH'
                %((z.imag + f_imag)/ self.omega * 10 ** 6 * self.z0.real))

            else :
                "Z is downside 1S admittance circle."

                # reactance of group 1 = 1/jwC
                self.s_r_1_t.set_text('____ s_Capacitor(C) : %.3f nF'
                %(1/(f_imag - z.imag)/ self.omega * 10 ** 9 / self.z0.real))
                # reactance of group 2 = 1/jwC
                self.s_r_3_t.set_text('____ s_Capacitor(C) : %.3f nF'
                %(1/(- f_imag - z.imag)/ self.omega * 10 ** 9 / self.z0.real))

        else :
            "Need to switch chart, first component is connected in parallel."

            # Equivalent admittance
            a = 1 / z

            # First susceptance with positive imaginary.
            f_imag = np.sqrt(a.real - a.real ** 2)
            # Second reactance with positive imaginary.
            s_imag = (-1 / complex(a.real, f_imag)).imag

            #Draw matching path.
            self.second_order_matching_series(a, f_imag)
            g = self.z2gamma(z)
            self.s_order_switch_t.set_data([g.real, -g.real], [g.imag, -g.imag])

            # Update reactance and susceptance label.

            # L of group 1 reactance = 1/jwc
            self.s_r_2_t.set_text('----- s_Capacitor(C) : %.3f nF'
            %(1 / s_imag / self.omega * 10 ** 9 / self.z0.real))
            # C of group 2 reactance = jwL
            self.s_r_4_t.set_text('----- s_Inductor(L) : %.3f uH'
            %(s_imag / self.omega * 10 ** 6 * self.z0.real))

            # susceptance of group 1 = 1/jwL
            self.s_r_1_t.set_text('____ p_Inductor(L) : %.3f uH'
            %(1/(f_imag - a.imag)/ self.omega * 10 ** 6 * self.z0.real))
            # susceptance of group 2 = jwC
            self.s_r_3_t.set_text('____ p_Capacitor(C) : %.3f nF'
            %((a.imag + f_imag)/ self.omega * 10 ** 9 / self.z0.real))

    def second_order_matching_series(self, z, t_imag):
        "Draw 2nd order impedance(normalized) matching path when it's cascade."

        tmp_imag = t_imag

        # Point lists to draw first-components path matching lines.
        reactance_1_1_l = []
        reactance_1_2_l = []

        # Two impedances with conductance equal to 1 transform to gamma form.
        g_tmp_1 = self.z2gamma(complex(z.real, tmp_imag))
        g_tmp_2 = self.z2gamma(complex(z.real, -tmp_imag))

        z_tmp_real = complex(1 / complex(z.real, -tmp_imag)).real
        z_tmp_imag = complex(1 / complex(z.real, -tmp_imag)).imag

        # Point lists to draw second-components path matching lines.
        reactance_2_1_l = [
        self.z2gamma(complex(z_tmp_real, y))
        for y in np.logspace(np.log10(z_tmp_imag), -2, 50)]
        reactance_2_2_l = [
        self.z2gamma(complex(z_tmp_real, -y))
        for y in np.logspace(-2, np.log10(z_tmp_imag), 50)]

        if z.imag > 0:
            reactance_1_1_l = [
                self.z2gamma(complex(z.real, y))
                for y in np.logspace(np.log10(tmp_imag),
                np.log10(z.imag), 100)]
            reactance_1_2_l = [
                self.z2gamma(complex(z.real, y))
                for y in np.logspace(np.log10(z.imag), -2, 50)] + [
                self.z2gamma(complex(z.real, -y))
                for y in np.logspace(-2, np.log10(tmp_imag), 50)]

        else :
            reactance_1_1_l = [
                self.z2gamma(complex(z.real, y))
                for y in np.logspace(np.log10(tmp_imag), -2, 50)] + [
                self.z2gamma(complex(z.real, -y))
                for y in np.logspace(-2, np.log10(-z.imag), 50)]
            reactance_1_2_l = [
                self.z2gamma(complex(z.real, -y))
                for y in np.logspace(np.log10(-z.imag),
                np.log10(tmp_imag), 100)]

        self.s_order_reactance_1_1.set_data([z.real for z in reactance_1_1_l],
            [z.imag for z in reactance_1_1_l])
        self.s_order_reactance_1_2.set_data([z.real for z in reactance_1_2_l],
            [z.imag for z in reactance_1_2_l])
        self.s_order_reactance_2_1.set_data([z.real for z in reactance_2_1_l],
            [z.imag for z in reactance_2_1_l])
        self.s_order_reactance_2_2.set_data([z.real for z in reactance_2_2_l],
            [z.imag for z in reactance_2_2_l])
        self.s_order_switch_1.set_data([g_tmp_1.real, -g_tmp_1.real],
            [g_tmp_1.imag, -g_tmp_1.imag])
        self.s_order_switch_2.set_data([g_tmp_2.real, -g_tmp_2.real],
            [g_tmp_2.imag, -g_tmp_2.imag])

    def show(self):
        """
        Shows the plot. The plot can't be updated after it has been closed.
        """
        plt.figure(self.fig.number)
        plt.show()

    def save(self, filename):
        """
        Saves the plot to filename. The extension defines the filetype.
        """
        self.fig.savefig(filename)

    def drawLine(self, spots, style, width):
        """
        Draws a list of spots and connects them by lines.
        """
        plt.figure(self.fig.number)
        xlst = [self.z2gamma(z).real for z in spots]
        ylst = [self.z2gamma(z).imag for z in spots]
        plt.plot(xlst, ylst, style, linewidth = width)
        plt.draw()

    def drawXCircle(self, x, style = 'k-', width = 1, p_num = 200):
        """
        Draws a circle with constant real part.
        """
        zlst = [x]+[complex(x,  y) for y in np.logspace(-3, 3, p_num)]
        self.drawLine(zlst, style, width)
        zlst = [x]+[complex(x, -y) for y in np.logspace(-3, 3, p_num)]
        self.drawLine(zlst, style, width)

    def drawYCircle(self, y, style = 'k-', width = 1, p_num = 200):
        """
        Draws a circle with constant imaginary part.
        """
        zlst = [complex(0,y)]+[complex(x, y) for x in np.logspace(-3, 3, p_num)]
        self.drawLine(zlst, style, width)

    def circleLable(self, txt, x, y, x_offset = 0, y_offset = 0):
        "Real or imaginary part value of a circle."

        c = self.z2gamma(complex(x,y))
        plt.figure(self.fig.number)
        plt.text(c.real + x_offset, c.imag + y_offset, txt, fontsize = 20)

    def z2gamma(self, z):
        """
        Converts an impedance to a reflection coefficient.
        'z' is normalized impedance.
        """
        return complex(z - 1)/(z + 1)

    def gamma2z(self, gamma):
        """
        Converts a reflection coefficient to impedance.
        'z' is normalized impedance.
        """
        return complex(1 + gamma)/(1 - gamma)


if __name__ == '__main__':
    s = Smith_co()
    s.show()
