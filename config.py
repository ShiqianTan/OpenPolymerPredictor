# ============= config.py =============
"""
Configuration settings for the Flask application
"""
import os

class Config:
    """Base configuration"""
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'dev-secret-key-change-in-production'
    DEBUG = False
    HOST = '0.0.0.0'
    PORT = 5000
    
    # Model configuration
    N_ESTIMATORS = 50
    RANDOM_STATE = 42
    FINGERPRINT_BITS = 2048
    MORGAN_RADIUS = 2

class DevelopmentConfig(Config):
    """Development configuration"""
    DEBUG = True

class ProductionConfig(Config):
    """Production configuration"""
    DEBUG = False
    SECRET_KEY = os.environ.get('SECRET_KEY')

class TestingConfig(Config):
    """Testing configuration"""
    TESTING = True
    DEBUG = True

# Default config
DEBUG = True
HOST = '0.0.0.0'
PORT = 5000
