
install ( FILES __init__.py IntegralCompoundPoisson.py IntegralUserDefinedFactory.py IntegralUserDefined.py Polynomial.py
          DESTINATION ${OTICP_PYTHON_MODULE_PATH}/${PACKAGE_NAME}
          PERMISSIONS OWNER_READ GROUP_READ WORLD_READ
        )
